/* $Id: test_ncbi_ftp_connector.c,v 1.24 2010/04/01 14:15:35 kazimird Exp $
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
 *   Test case for FTP-based CONNECTOR
 *
 */

#include <connect/ncbi_connection.h>
#include <connect/ncbi_connutil.h>
#include <connect/ncbi_ftp_connector.h>
#include "../ncbi_ansi_ext.h"
#include "../ncbi_priv.h"               /* CORE logging facilities */
#include <stdlib.h>
#ifdef HAVE_GETTIMEOFDAY
#  include <sys/time.h>
#endif /*HAVE_GETTIMEOFDAY*/
#include <time.h>
/* This header must go last */
#include "test_assert.h"

#define TEST_HOST "ftp.ncbi.nlm.nih.gov"
#define TEST_PORT 0
#define TEST_USER "ftp"
#define TEST_PASS "none"
#define TEST_PATH ((char*) 0)


static double s_GetTime(void)
{
#ifdef HAVE_GETTIMEOFDAY
    struct timeval t;
    return gettimeofday(&t, 0) == 0 ? t.tv_sec + t.tv_usec / 1000000.0 : 0.0;
#else
    time_t t = time(0);
    return (double)((unsigned long) t);
#endif /*HAVE_GETTIMEOFDAY*/
}


int main(int argc, char* argv[])
{
    static const char kChdir[] = "CWD /toolbox/ncbi_tools\n";
    static const char kFile[] = "RETR CURRENT/ncbi.tar.gz";
    int/*bool*/   aborting = 0, first;
    TFCDC_Flags   flag = 0;
    SConnNetInfo* net_info;
    char          buf[1024];
    CONNECTOR     connector;
    FILE*         data_file;
    size_t        size, n;
    EIO_Status    status;
    double        elapse;
    CONN          conn;

    g_NCBI_ConnectRandomSeed = (int) time(0) ^ NCBI_CONNECT_SRAND_ADDEND;
    srand(g_NCBI_ConnectRandomSeed);

    /* Log and data-log streams */
    CORE_SetLOGFormatFlags(fLOG_None          | fLOG_Level   |
                           fLOG_OmitNoteLevel | fLOG_DateTime);
    CORE_SetLOGFILE(stderr, 0/*false*/);
    data_file = fopen("test_ncbi_ftp_connector.dat", "wb");
    assert(data_file);

    assert((net_info = ConnNetInfo_Create(0)) != 0);
    if (net_info->debug_printout == eDebugPrintout_Some)
        flag |= fFCDC_LogControl;
    else if (net_info->debug_printout == eDebugPrintout_Data) {
        char val[32];
        ConnNetInfo_GetValue(0, REG_CONN_DEBUG_PRINTOUT, val, sizeof(val),
                             DEF_CONN_DEBUG_PRINTOUT);
        flag |= strcasecmp(val, "all") == 0 ? fFCDC_LogAll : fFCDC_LogData;
    }

    if (TEST_PORT) {
        sprintf(buf, ":%hu", TEST_PORT);
    } else {
        *buf = 0;
    }
    CORE_LOGF(eLOG_Note, ("Connecting to ftp://%s:%s@%s%s/",
                          TEST_USER, TEST_PASS, TEST_HOST, buf));
    /* Run the tests */
    connector = FTP_CreateDownloadConnector(TEST_HOST, TEST_PORT,
                                            TEST_USER, TEST_PASS,
                                            TEST_PATH, flag);

    if (CONN_Create(connector, &conn) != eIO_Success)
        CORE_LOG(eLOG_Fatal, "Cannot create FTP download connection");

    assert(CONN_SetTimeout(conn, eIO_Open,      net_info->timeout)
           == eIO_Success);
    assert(CONN_SetTimeout(conn, eIO_ReadWrite, net_info->timeout)
           == eIO_Success);
    assert(CONN_SetTimeout(conn, eIO_Close,     net_info->timeout)
           == eIO_Success);

    if (CONN_Read(conn, buf, sizeof(buf), &n, eIO_ReadPlain) != eIO_Closed)
        CORE_LOG(eLOG_Fatal, "Test failed in empty READ");

    if (CONN_Write(conn, "aaa", 3, &n, eIO_WritePlain) != eIO_Success)
        CORE_LOG(eLOG_Fatal, "Cannot write FTP command");

    if (CONN_Wait(conn, eIO_Read, net_info->timeout) != eIO_Unknown)
        CORE_LOG(eLOG_Fatal, "Test failed in waiting on READ");
    CORE_LOG(eLOG_Note, "Unrecognized command was correctly rejected");

    if (CONN_Write(conn, "LIST\nSIZE", 9, &n, eIO_WritePlain) != eIO_Unknown)
        CORE_LOG(eLOG_Fatal, "Test failed to reject multiple commands");
    CORE_LOG(eLOG_Note, "Multiple commands were correctly rejected");

    if (CONN_Write(conn, "LIST", 4, &n, eIO_WritePlain) != eIO_Success)
        CORE_LOG(eLOG_Fatal, "Cannot write LIST command");

    CORE_LOG(eLOG_Note, "LIST command output:");
    first = 1/*true*/;
    do {
        status = CONN_Read(conn, buf, sizeof(buf), &n, eIO_ReadPlain);
        if (n != 0) {
            printf("%.*s", (int) n, buf);
            first = 0/*false*/;
            fflush(stdout);
        }
    } while (status == eIO_Success);
    if (first  ||  status != eIO_Closed) {
        printf("<%s>\n", status != eIO_Success ? IO_StatusStr(status) : "EOF");
        fflush(stdout);
    }

    if (CONN_Write(conn, "NLST\r\n", 6, &n, eIO_WritePlain) != eIO_Success)
        CORE_LOG(eLOG_Fatal, "Cannot write NLST command");

    CORE_LOG(eLOG_Note, "NLST command output:");
    first = 1/*true*/;
    do {
        status = CONN_Read(conn, buf, sizeof(buf), &n, eIO_ReadPlain);
        if (n != 0) {
            printf("%.*s", (int) n, buf);
            first = 0/*false*/;
            fflush(stdout);
        } else {
            assert(status != eIO_Success);
        }
    } while (status == eIO_Success);
    if (first  ||  status != eIO_Closed) {
        printf("<%s>\n", status != eIO_Success ? IO_StatusStr(status) : "EOF");
        fflush(stdout);
    }

    if (CONN_Write(conn, kChdir, sizeof(kChdir) - 1, &n, eIO_WritePlain)
        != eIO_Success) {
        CORE_LOGF(eLOG_Fatal, ("Cannot execute %.*s",
                               (int) sizeof(kChdir) - 2, kChdir));
    }

    if (CONN_Write(conn, kFile, sizeof(kFile) - 1, &n, eIO_WritePersist)
        != eIO_Success) {
        CORE_LOGF(eLOG_Fatal, ("Cannot write %s", kFile));
    }

    size = 0;
    elapse = s_GetTime();
    do {
        status = CONN_Read(conn, buf, sizeof(buf), &n, eIO_ReadPlain);
        if (n != 0) {
            fwrite(buf, n, 1, data_file);
            fflush(data_file);
            size += n;
            rand();
            if (argc > 1  &&  rand() % 100 == 0) {
                aborting = 1;
                break;
            }
        } else {
            assert(status != eIO_Success);
            if (status != eIO_Closed  ||  !size)
                CORE_LOGF(eLOG_Error, ("Read error: %s",IO_StatusStr(status)));
        }
    } while (status == eIO_Success);
    elapse = s_GetTime() - elapse;

    if (!aborting  ||  (rand() & 1) == 0) {
        if (CONN_Write(conn, "NLST blah*", 10, &n, eIO_WritePlain)
            != eIO_Success) {
            CORE_LOG(eLOG_Fatal, "Cannot write garbled NLST command");
        }

        CORE_LOG(eLOG_Note, "Garbled NLST command output (should be empty):");
        first = 1/*true*/;
        do {
            status = CONN_Read(conn, buf, sizeof(buf), &n, eIO_ReadPlain);
            if (n != 0) {
                printf("%.*s", (int) n, buf);
                first = 0/*false*/;
                fflush(stdout);
            }
        } while (status == eIO_Success);
        if (first) {
            printf("<EOF>\n");
            fflush(stdout);
        }
    }

    if (CONN_Close(conn) != eIO_Success) {
        CORE_LOGF(eLOG_Fatal, ("Error %s FTP connection",
                               aborting ? "aborting" : "closing"));
    }

    /* Cleanup and exit */
    fclose(data_file);
    if (!aborting) {
        CORE_LOGF(size ? eLOG_Note : eLOG_Fatal,
                  ("%lu byte(s) downloaded in %.2f second(s) @ %.2fKB/s",
                   (unsigned long) size, elapse,
                   (unsigned long) size / (1024 * (elapse ? elapse : 1.0))));
    } else
        remove("test_ncbi_ftp_connector.dat");
    ConnNetInfo_Destroy(net_info);

    CORE_LOG(eLOG_Note, "TEST completed successfully");
    CORE_SetLOG(0);
    return 0;
}
