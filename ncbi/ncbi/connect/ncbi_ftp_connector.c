/* $Id: ncbi_ftp_connector.c,v 1.33 2010/03/01 13:54:31 kazimird Exp $
 * ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * ===========================================================================
 *
 * Author:  Anton Lavrentiev
 *
 * File Description:
 *   FTP CONNECTOR
 *   See also:  RFCs 959 (STD 9), 1634 (FYI 24),
 *   and IETF 9-2002 "Extensions to FTP".
 *
 *   See <connect/ncbi_connector.h> for the detailed specification of
 *   the connector's methods and structures.
 *
 */

#include "ncbi_ansi_ext.h"
#include "ncbi_assert.h"
#include "ncbi_priv.h"
#include <connect/ncbi_ftp_connector.h>
#include <connect/ncbi_socket.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>

#define NCBI_USE_ERRCODE_X   Connect_FTP


/***********************************************************************
 *  INTERNAL -- Auxiliary types and static functions
 ***********************************************************************/

typedef enum {
    fFtpFeature_MDTM = 1,
    fFtpFeature_SIZE = 2
} EFTP_Feature;
typedef unsigned int TFTP_Features; /* bitwise OR of EFtpFeature's */


/* All internal data necessary to perform I/O
 */
typedef struct {
    const char*    host;
    unsigned short port;
    const char*    user;
    const char*    pass;
    const char*    path;
    const char*    name;
    TFTP_Features  feat;
    TFCDC_Flags    flag;
    SOCK           cntl;     /* control connection */
    SOCK           data;     /* data    connection */
    BUF            wbuf;     /* write buffer       */
    EIO_Status     r_status; /* pertains to data    connection IO */
    EIO_Status     w_status; /* pertains to control connection IO */
} SFTPConnector;


static const STimeout kFTPFailsafeTimeout = {10, 0};


static EIO_Status s_ReadReply(SOCK sock, int* code,
                              char* line, size_t maxlinelen)
{
    int/*bool*/ first = 1/*true*/;
    for (;;) {
        int c, m;
        size_t n;
        char buf[1024];
        EIO_Status status = SOCK_ReadLine(sock, buf, sizeof(buf), &n);
        /* All FTP replies are at least '\n'-terminated, no ending with EOF */
        if (status != eIO_Success)
            return status;
        if (n == sizeof(buf))
            return eIO_Unknown/*line too long*/;
        if (first  ||  isdigit((unsigned char) *buf)) {
            if (sscanf(buf, "%d%n", &c, &m) < 1)
                return eIO_Unknown;
        } else
            c = 0;
        if (first) {
            if (m != 3  ||  code == 0)
                return eIO_Unknown;
            if (line)
                strncpy0(line, &buf[m + 1], maxlinelen);
            *code = c;
            if (buf[m] != '-') {
                if (buf[m] == ' ')
                    break;
                return eIO_Unknown;
            }
            first = 0/*false*/;
        } else if (c == *code  &&  m == 3  &&  buf[m] == ' ')
            break;
    }
    return eIO_Success;
}


static EIO_Status s_FTPCloseData(SFTPConnector* xxx, int/*bool*/ abort)
{
    EIO_Status status;
    assert(xxx->data);
    if (xxx->flag & fFCDC_LogControl)
        SOCK_SetDataLogging(xxx->data, eOn);
    if (!abort) {
        static const STimeout kInstant = {0, 0};
        SOCK_SetTimeout(xxx->data, eIO_Close, &kInstant);
        status = eIO_Success;
    } else
        status = SOCK_Abort(xxx->data);
    SOCK_Close(xxx->data);
    xxx->data = 0;
    return status;
}


static EIO_Status s_FTPReply(SFTPConnector* xxx, int* code,
                             char* line, size_t maxlinelen)
{
    int c = 0;
    EIO_Status status = eIO_Closed;
    if (xxx->cntl) {
        status = s_ReadReply(xxx->cntl, &c, line, maxlinelen);
        if (status == eIO_Success  &&  c == 421)
            status =  eIO_Closed;
        if (status == eIO_Closed  ||  (status == eIO_Success  &&  c == 221)) {
            if (xxx->data) {
                s_FTPCloseData(xxx,
                               status == eIO_Closed ? 1/*abort*/ : 0/*close*/);
            }
            SOCK_Close(xxx->cntl);
            xxx->cntl = 0;
        }
    }
    if (code)
        *code = c;
    return status;
}


static EIO_Status s_FTPDrainReply(SFTPConnector* xxx, int* code, int cXX)
{
    EIO_Status status;
    int        c;
    while ((status = s_FTPReply(xxx, &c, 0, 0)) == eIO_Success
           &&  (!cXX  ||  c/100 != cXX)) {
        *code = c;
    }
    return status;
}


static EIO_Status s_FTPCommand(SFTPConnector* xxx,
                               const char* cmd, const char* arg)
{
    size_t cmdlen  = strlen(cmd);
    size_t arglen  = arg ? strlen(arg) : 0;
    size_t linelen = cmdlen + (arglen ? 1/* */ + arglen : 0) + 2/*\r\n*/;
    char*  line    = (char*) malloc(linelen + 1/*\0*/);
    EIO_Status status;

    assert(xxx->cntl);
    if (line) {
        memcpy(line, cmd, cmdlen);
        if (arglen) {
            line[cmdlen++] = ' ';
            memcpy(line + cmdlen, arg, arglen);
            cmdlen += arglen;
        }
        line[cmdlen++] = '\r';
        line[cmdlen++] = '\n';
        line[cmdlen]   = '\0';
        status = SOCK_Write(xxx->cntl, line, linelen, 0, eIO_WritePersist);
        free(line);
    } else
        status = eIO_Unknown;
    return status;
}


static EIO_Status s_FTPLogin(SFTPConnector* xxx, const STimeout* timeout)
{
    EIO_Status status;
    int code;

    assert(xxx->cntl  &&  xxx->user  &&  xxx->pass);
    SOCK_SetTimeout(xxx->cntl, eIO_ReadWrite, timeout);
    status = s_FTPReply(xxx, &code, 0, 0);
    if (status != eIO_Success)
        return status;
    if (code != 220  ||  !*xxx->user)
        return eIO_Unknown;
    status = s_FTPCommand(xxx, "USER", xxx->user);
    if (status != eIO_Success)
        return status;
    status = s_FTPReply(xxx, &code, 0, 0);
    if (status != eIO_Success)
        return status;
    if (code == 230)
        return eIO_Success;
    if (code != 331  ||  !*xxx->pass)
        return eIO_Unknown;
    status = s_FTPCommand(xxx, "PASS", xxx->pass);
    if (status != eIO_Success)
        return status;
    status = s_FTPReply(xxx, &code, 0, 0);
    if (status != eIO_Success)
        return status;
    return code == 230 ? eIO_Success : eIO_Unknown;
}


static EIO_Status s_FTPChdir(SFTPConnector* xxx, const char* cmd)
{
    if (cmd  ||  xxx->path) {
        int code;
        EIO_Status status = s_FTPCommand(xxx,
                                         cmd ? cmd : "CWD",
                                         cmd ? 0   : xxx->path);
        if (status != eIO_Success)
            return status;
        status = s_FTPReply(xxx, &code, 0, 0);
        if (status != eIO_Success)
            return status;
        if (code != 250)
            return eIO_Unknown;
    }
    return eIO_Success;
}


static EIO_Status s_FTPBinary(SFTPConnector* xxx)
{
    int code;
    EIO_Status status = s_FTPCommand(xxx, "TYPE", "I");
    if (status != eIO_Success)
        return status;
    status = s_FTPReply(xxx, &code, 0, 0);
    if (status != eIO_Success)
        return status;
    return code == 200 ? eIO_Success : eIO_Unknown;
}


static EIO_Status s_FTPAbort(SFTPConnector*  xxx,
                             const STimeout* timeout,
                             int/*bool*/     quit)
{
    EIO_Status status = eIO_Success;
    int        code;
    size_t     n;

    if (!xxx->data)
        return status;
    if (quit  ||  !xxx->cntl)
        return s_FTPCloseData(xxx, 1/*abort*/);
    if (!timeout)
        timeout = &kFTPFailsafeTimeout;
    SOCK_SetTimeout(xxx->cntl, eIO_ReadWrite, timeout);
    SOCK_SetTimeout(xxx->data, eIO_ReadWrite, timeout);
    if (/* Send TELNET IP (Interrupt Process) command */
        (status = SOCK_Write(xxx->cntl, "\377\364", 2, &n, eIO_WritePersist))
        != eIO_Success  ||  n != 2                                         ||
        /* Send TELNET DM (Data Mark) command to complete SYNCH, RFC 854 */
        (status = SOCK_Write(xxx->cntl, "\377\362", 2, &n, eIO_WriteOutOfBand))
        != eIO_Success  ||  n != 2                                         ||
        (status = s_FTPCommand(xxx, "ABOR", 0)) != eIO_Success) {
        s_FTPCloseData(xxx, 1/*abort*/);
        return status == eIO_Success ? eIO_Unknown : status;
    }
    while (SOCK_Read(xxx->data, 0, 1024*1024/*drain up*/, 0, eIO_ReadPlain)
           == eIO_Success) {
        continue;
    }
    n = status == eIO_Timeout ? 1 : 0;
    s_FTPCloseData(xxx, SOCK_Status(xxx->data, eIO_Read) == eIO_Closed
                   ? 0/*close*/ : 1/*abort*/);
    if ((status = s_FTPDrainReply(xxx, &code, 2/*2xx*/)) == eIO_Success
        /* Microsoft FTP is known to return 225 (instead of 226) */
        &&  code != 225  &&  code != 226) {
        status = eIO_Unknown;
    }
    if (n  ||  status == eIO_Timeout) {
        CORE_LOG_X(1, eLOG_Warning,
                   "[FTP]  Timed out in data connection abort");
    }
    return status;
}


static EIO_Status s_FTPPasv(SFTPConnector* xxx)
{
    static const STimeout kInstant = {0, 0};
    unsigned int   host, i;
    unsigned short port;
    EIO_Status status;
    int  code, o[6];
    char buf[128];

    status = s_FTPCommand(xxx, "PASV", 0);
    if (status != eIO_Success)
        return status;
    status = s_FTPReply(xxx, &code, buf, sizeof(buf) - 1);
    if (status != eIO_Success  ||  code != 227)
        return eIO_Unknown;
    buf[sizeof(buf) - 1] = '\0';
    for (;;) {
        char* c;
        /* RFC 1123 4.1.2.6 says that ()'s in PASV reply must not be assumed */
        for (c = buf; *c; c++) {
            if (isdigit((unsigned char)(*c)))
                break;
        }
        if (!*c)
            return eIO_Unknown;
        if (sscanf(c, "%d,%d,%d,%d,%d,%d%n",
                   &o[0], &o[1], &o[2], &o[3], &o[4], &o[5], &code) >= 6) {
            break;
        }
        memmove(buf, c + code, strlen(c + code) + 1);
    }
    for (i = 0;  i < (unsigned int)(sizeof(o)/sizeof(o[0]));  i++) {
        if (o[i] < 0  ||  o[i] > 255)
            return eIO_Unknown;
    }
    i = (((((o[0] << 8) | o[1]) << 8) | o[2]) << 8) | o[3];
    host = SOCK_HostToNetLong(i);
    i = (o[4] << 8) | o[5];
    port = (unsigned short) i;
    if (SOCK_ntoa(host, buf, sizeof(buf)) != 0)
        return eIO_Unknown;
    status = SOCK_CreateEx(buf, port, &kInstant, &xxx->data, 0, 0,
                           xxx->flag & fFCDC_LogControl
                           ? fSOCK_LogOn : fSOCK_LogDefault);
    if (status != eIO_Success) {
        CORE_LOGF_X(2, eLOG_Error,
                    ("[FTP]  Cannot open data connection to %s:%hu: %s",
                     buf, port, IO_StatusStr(status)));
        s_FTPAbort(xxx, 0, 0/*!quit*/);
        assert(!xxx->data);
        return eIO_Unknown;
    }
    if (!(xxx->flag & fFCDC_LogData))
        SOCK_SetDataLogging(xxx->data, eDefault);
    return eIO_Success;
}


static EIO_Status s_FTPRetrieve(SFTPConnector* xxx,
                                const char*    cmd)
{
    int code;
    EIO_Status status = s_FTPPasv(xxx);
    if (status != eIO_Success)
        return status;
    status = s_FTPCommand(xxx, cmd, 0);
    if (status != eIO_Success)
        return status;
    status = s_FTPReply(xxx, &code, 0, 0);
    if (status != eIO_Success)
        return status;
    if (code == 150)
        return eIO_Success;
    if (code == 450  &&  (strncasecmp(cmd, "NLST", 4) == 0  ||
                          strncasecmp(cmd, "LIST", 4) == 0)) {
        /* server usually drops data connection on 450: no files ...*/
        if (xxx->data)
            s_FTPCloseData(xxx, 1/*abort*/);
        /* with no data connection open, user gets eIO_Closed on read */
        return eIO_Success;
    }
    s_FTPAbort(xxx, 0, 0/*!quit*/);
    assert(!xxx->data);
    return eIO_Unknown;
}


static EIO_Status s_FTPExecute(SFTPConnector* xxx, const STimeout* timeout)
{
    EIO_Status status;
    size_t     size;
    char*      s;

    status = s_FTPAbort(xxx, timeout, 0/*!quit*/);
    assert(!xxx->data);
    if (status != eIO_Success)
        return status;
    assert(xxx->cntl);
    verify((size = BUF_Size(xxx->wbuf)) > 0);
    if (!(s = (char*) malloc(size + 1)))
        return eIO_Unknown;
    if (BUF_Read(xxx->wbuf, s, size) == size) {
        const char* c;
        assert(!memchr(s, '\n', size));
        if (s[size - 1] == '\r')
            --size;
        s[size] = '\0';
        if (!(c = (const char*) memchr(s, ' ', size)))
            c = s + size;
        else
            size = (size_t)(c - s);
        if (size == 3  ||  size == 4) {
            SOCK_SetTimeout(xxx->cntl, eIO_ReadWrite, timeout);
            if        (size == 3  &&  strncasecmp(s, "CWD",  3) == 0) {
                status = s_FTPChdir(xxx, s);
            } else if (size == 4  && (strncasecmp(s, "LIST", 4) == 0  ||
                                      strncasecmp(s, "NLST", 4) == 0  ||
                                      strncasecmp(s, "RETR", 4) == 0)) {
                status = s_FTPRetrieve(xxx, s);
            } else if (size == 4  &&  strncasecmp(s, "REST", 4) == 0) {
                status = s_FTPCommand(xxx, s, 0);
                if (status == eIO_Success) {
                    int code;
                    status = s_FTPReply(xxx, &code, 0, 0);
                    if (status == eIO_Success  &&  code != 350)
                        status = eIO_Unknown;
                }
            } else
                status = eIO_Unknown;
        } else
            status = eIO_Unknown;
    } else
        status = eIO_Unknown;
    free(s);
    return status;
}


/***********************************************************************
 *  INTERNAL -- "s_VT_*" functions for the "virt. table" of connector methods
 ***********************************************************************/

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
    static const char* s_VT_GetType (CONNECTOR       connector);
    static EIO_Status  s_VT_Open    (CONNECTOR       connector,
                                     const STimeout* timeout);
    static EIO_Status  s_VT_Wait    (CONNECTOR       connector,
                                     EIO_Event       event,
                                     const STimeout* timeout);
    static EIO_Status  s_VT_Write   (CONNECTOR       connector,
                                     const void*     buf,
                                     size_t          size,
                                     size_t*         n_written,
                                     const STimeout* timeout);
    static EIO_Status  s_VT_Read    (CONNECTOR       connector,
                                     void*           buf,
                                     size_t          size,
                                     size_t*         n_read,
                                     const STimeout* timeout);
    static EIO_Status  s_VT_Flush   (CONNECTOR       connector,
                                     const STimeout* timeout);
    static EIO_Status  s_VT_Status  (CONNECTOR       connector,
                                     EIO_Event       dir);
    static EIO_Status  s_VT_Close   (CONNECTOR       connector,
                                     const STimeout* timeout);
    static void        s_Setup      (SMetaConnector* meta,
                                     CONNECTOR       connector);
    static void        s_Destroy    (CONNECTOR       connector);
#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */


/*ARGSUSED*/
static const char* s_VT_GetType
(CONNECTOR connector)
{
    return "FTP";
}


static EIO_Status s_VT_Open
(CONNECTOR       connector,
 const STimeout* timeout)
{
    SFTPConnector* xxx = (SFTPConnector*) connector->handle;
    EIO_Status status;

    assert(!xxx->data  &&  !xxx->cntl);

    status = SOCK_CreateEx(xxx->host, xxx->port, timeout, &xxx->cntl, 0, 0,
                           xxx->flag & fFCDC_LogControl
                           ? fSOCK_LogOn : fSOCK_LogDefault);
    if (status == eIO_Success) {
        SOCK_DisableOSSendDelay(xxx->cntl, 1/*yes,disable*/);
        status =  s_FTPLogin(xxx, timeout);
    }
    if (status == eIO_Success)
        status =  s_FTPChdir(xxx, 0);
    if (status == eIO_Success)
        status =  s_FTPBinary(xxx);
    if (status != eIO_Success) {
        if (xxx->cntl) {
            SOCK_Close(xxx->cntl);
            xxx->cntl = 0;
        }
    }
    xxx->r_status = status;
    xxx->w_status = status;
    return status;
}
 

static EIO_Status s_VT_Wait
(CONNECTOR       connector,
 EIO_Event       event,
 const STimeout* timeout)
{
    SFTPConnector* xxx = (SFTPConnector*) connector->handle;
    assert(event == eIO_Read  ||  event == eIO_Write);

    if (!xxx->cntl)
        return eIO_Closed;
    if (event & eIO_Write)
        return eIO_Success;
    if (!xxx->data) {
        EIO_Status status;
        if (!BUF_Size(xxx->wbuf))
            return eIO_Closed;
        if ((status = s_FTPExecute(xxx, timeout)) != eIO_Success)
            return status;
        if (!xxx->data)
            return eIO_Closed;
    }
    return SOCK_Wait(xxx->data, eIO_Read, timeout);
}


static EIO_Status s_VT_Write
(CONNECTOR       connector,
 const void*     buf,
 size_t          size,
 size_t*         n_written,
 const STimeout* timeout)
{
    SFTPConnector* xxx = (SFTPConnector*) connector->handle;
    EIO_Status status;
    const char* c;

    if (!xxx->cntl)
        return eIO_Closed;
    if (!size)
        return eIO_Success;
    if (((c = (const char*) memchr((const char*) buf, '\n', size)) != 0
         &&  c < (const char*) buf + size - 1)
        ||  !BUF_Write(&xxx->wbuf, buf, size - (c ? 1 : 0))) {
        return eIO_Unknown;
    }

    status = c ? s_FTPExecute(xxx, timeout) : eIO_Success;
    if (status == eIO_Success)
        *n_written = size;
    if (status != eIO_Timeout)
        xxx->w_status = status;
    return status;
}


static EIO_Status s_VT_Flush
(CONNECTOR       connector,
 const STimeout* timeout)
{
    SFTPConnector* xxx = (SFTPConnector*) connector->handle;
    return !xxx->cntl ? eIO_Closed
        : BUF_Size(xxx->wbuf) ? s_FTPExecute(xxx, timeout) : eIO_Success;
}


static EIO_Status s_VT_Read
(CONNECTOR       connector,
 void*           buf,
 size_t          size,
 size_t*         n_read,
 const STimeout* timeout)
{
    SFTPConnector* xxx = (SFTPConnector*) connector->handle;
    EIO_Status status = eIO_Closed;

    if (xxx->cntl  &&  xxx->data) {
        SOCK_SetTimeout(xxx->data, eIO_Read, timeout);
        status = SOCK_Read(xxx->data, buf, size, n_read, eIO_ReadPlain);
        if (status == eIO_Closed) {
            int code;
            s_FTPCloseData(xxx, 0/*close*/);
            SOCK_SetTimeout(xxx->cntl, eIO_Read, timeout);
            if (s_FTPReply(xxx, &code, 0, 0) != eIO_Success
                ||  (code != 225  &&  code != 226)) {
                status = eIO_Unknown;
            }
        }
        if (status != eIO_Timeout)
            xxx->r_status = status;
    }
    return status;
}


static EIO_Status s_VT_Status
(CONNECTOR connector,
 EIO_Event dir)
{
    SFTPConnector* xxx = (SFTPConnector*) connector->handle;

    switch (dir) {
    case eIO_Read:
        return xxx->data ? xxx->r_status : eIO_Closed;
    case eIO_Write:
        return xxx->cntl ? xxx->w_status : eIO_Closed;
    default:
        assert(0); /* should never happen as checked by connection */
        break;
    }
    return eIO_InvalidArg;
}


static EIO_Status s_VT_Close
(CONNECTOR       connector,
 const STimeout* timeout)
{
    SFTPConnector* xxx = (SFTPConnector*) connector->handle;
    EIO_Status status;

    status = s_FTPAbort(xxx, timeout, 1/*quit*/);
    assert(!xxx->data);
    if (status == eIO_Success) {
        if (xxx->cntl) {
            int code;
            if (!timeout)
                timeout = &kFTPFailsafeTimeout;
            SOCK_SetTimeout(xxx->cntl, eIO_ReadWrite, timeout);
            status = s_FTPCommand(xxx, "QUIT", 0);
            if (status == eIO_Success)
                status = s_FTPDrainReply(xxx, &code, 0);
            if (status != eIO_Closed  ||  code != 221)
                status = eIO_Unknown;
        }
    }
    if (xxx->cntl) {
        assert(status != eIO_Success);
        if (status == eIO_Timeout)
            SOCK_Abort(xxx->cntl);
        status = SOCK_Close(xxx->cntl);
        xxx->cntl = 0;
    }
    return status;
}


static void s_Setup
(SMetaConnector* meta,
 CONNECTOR       connector)
{
    /* initialize virtual table */
    CONN_SET_METHOD(meta, get_type, s_VT_GetType, connector);
    CONN_SET_METHOD(meta, open,     s_VT_Open,    connector);
    CONN_SET_METHOD(meta, wait,     s_VT_Wait,    connector);
    CONN_SET_METHOD(meta, write,    s_VT_Write,   connector);
    CONN_SET_METHOD(meta, flush,    s_VT_Flush,   connector);
    CONN_SET_METHOD(meta, read,     s_VT_Read,    connector);
    CONN_SET_METHOD(meta, status,   s_VT_Status,  connector);
    CONN_SET_METHOD(meta, close,    s_VT_Close,   connector);
    meta->default_timeout = kInfiniteTimeout;
}


static void s_Destroy
(CONNECTOR connector)
{
    SFTPConnector* xxx = (SFTPConnector*) connector->handle;
    connector->handle = 0;

    if (xxx->host) {
        free((void*) xxx->host);
        xxx->host = 0;
    }
    if (xxx->user) {
        free((void*) xxx->user);
        xxx->user = 0;
    }
    if (xxx->pass) {
        free((void*) xxx->pass);
        xxx->pass = 0;
    }
    if (xxx->path) {
        free((void*) xxx->path);
        xxx->path = 0;
    }
    if (xxx->name) {
        free((void*) xxx->name);
        xxx->name = 0;
    }
    BUF_Destroy(xxx->wbuf);
    xxx->wbuf = 0;
    free(xxx);
    free(connector);
}


/***********************************************************************
 *  EXTERNAL -- the connector's "constructors"
 ***********************************************************************/

extern CONNECTOR FTP_CreateDownloadConnector(const char*    host,
                                             unsigned short port,
                                             const char*    user,
                                             const char*    pass,
                                             const char*    path,
                                             TFCDC_Flags    flag)
{
    CONNECTOR      ccc = (SConnector*)    malloc(sizeof(SConnector));
    SFTPConnector* xxx = (SFTPConnector*) malloc(sizeof(*xxx));

    assert(!(flag & ~fFCDC_LogAll));

    xxx->data     = 0;
    xxx->cntl     = 0;
    xxx->wbuf     = 0;
    xxx->host     = strdup(host);
    xxx->port     = port ? port : 21;
    xxx->user     = strdup(user ? user : "ftp");
    xxx->pass     = strdup(pass ? pass : "-none");
    xxx->path     = path  &&  *path ? strdup(path) : 0;
    xxx->name     = 0;
    xxx->flag     = flag;

    /* initialize connector data */
    ccc->handle   = xxx;
    ccc->next     = 0;
    ccc->meta     = 0;
    ccc->setup    = s_Setup;
    ccc->destroy  = s_Destroy;

    return ccc;
}
