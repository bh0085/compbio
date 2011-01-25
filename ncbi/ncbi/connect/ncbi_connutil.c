/* $Id: ncbi_connutil.c,v 6.146 2010/05/03 04:05:41 kazimird Exp $
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
 * Author:  Denis Vakatov, Anton Lavrentiev
 *
 * File Description:
 *   Auxiliary API, mostly CONN-, URL-, and MIME-related
 *   (see in "ncbi_connutil.h" for more details).
 *
 */

#include "ncbi_ansi_ext.h"
#include "ncbi_priv.h"
#include <connect/ncbi_connutil.h>
#include <ctype.h>
#include <errno.h>
#include <stdlib.h>

#define NCBI_USE_ERRCODE_X   Connect_Util


static char* x_StrcatCRLF(char* dst, const char* src)
{
    size_t dstlen = dst  &&  *dst ? strlen(dst) : 0;
    size_t srclen = src  &&  *src ? strlen(src) : 0;
    if (dstlen  ||  srclen) {
        size_t len;
        char*  temp;
        if (dstlen  &&  dst[dstlen - 1] == '\n') {
            if (--dstlen  &&  dst[dstlen - 1] == '\r')
                --dstlen;
        }
        if (srclen  &&  src[srclen - 1] == '\n') {
            if (--srclen  &&  src[srclen - 1] == '\r')
                --srclen;
        }
        len = (dstlen ? dstlen + 2 : 0) + (srclen ? srclen + 2 : 0) + 1;
        if (!(temp = (char*)(dst ? realloc(dst, len) : malloc (len))))
            return 0;
        dst = temp;
        if (dstlen) {
            temp += dstlen;
            memcpy(temp, "\r\n", 3);
            temp += 2;
        }
        if (srclen) {
            memcpy(temp, src, srclen);
            temp += srclen;
            memcpy(temp, "\r\n", 3);
        }
    }
    return dst;
}


extern const char* ConnNetInfo_GetValue(const char* service, const char* param,
                                        char* value, size_t value_size,
                                        const char* def_value)
{
    char        buf[128];
    const char* val;
    char*       s;

    if (!value  ||  value_size <= 0)
        return 0;
    *value = '\0';

    if (!param  ||  !*param)
        return 0;

    if (service  &&  *service) {
        /* Service-specific inquiry */
        char   temp[sizeof(buf)];
        size_t slen = strlen(service);
        size_t plen = strlen(param) + 1;
        size_t len  = slen + 1 + sizeof(DEF_CONN_REG_SECTION) + plen;
        if (len > sizeof(buf))
            return 0;

        /* First, environment search for 'service_CONN_param' */
        s = (char*) memcpy(buf, service, slen) + slen;
        *s++ = '_';
        memcpy(s, DEF_CONN_REG_SECTION, sizeof(DEF_CONN_REG_SECTION) - 1);
        s += sizeof(DEF_CONN_REG_SECTION) - 1;
        *s++ = '_';
        memcpy(s, param, plen);
        if ((val = getenv(strupr((char*) memcpy(temp, buf, len)))) != 0
            ||  (memcmp(temp, buf, len) != 0  &&  (val = getenv(buf)) != 0)) {
            return strncpy0(value, val, value_size - 1);
        }

        /* Next, search for 'CONN_param' in '[service]' registry section */
        buf[slen++] = '\0';
        s = buf + slen;
        CORE_REG_GET(buf, s, value, value_size, 0);
        if (*value)
            return value;
    } else {
        /* Common case. Form 'CONN_param' */
        size_t plen = strlen(param) + 1;
        if (sizeof(DEF_CONN_REG_SECTION) + plen > sizeof(buf))
            return 0;
        s = buf;
        memcpy(s, DEF_CONN_REG_SECTION, sizeof(DEF_CONN_REG_SECTION) - 1);
        s += sizeof(DEF_CONN_REG_SECTION) - 1;
        *s++ = '_';
        memcpy(s, param, plen);
        s = strupr(buf);
    }

    /* Environment search for 'CONN_param' */
    if ((val = getenv(s)) != 0)
        return strncpy0(value, val, value_size - 1);

    /* Last resort: Search for 'param' in default registry section */
    s += sizeof(DEF_CONN_REG_SECTION);
    CORE_REG_GET(DEF_CONN_REG_SECTION, s, value, value_size, def_value);
    return value;
}


int/*bool*/ ConnNetInfo_Boolean(const char* str)
{
    return str  &&  *str  &&  (strcmp    (str, "1")    == 0  ||
                               strcasecmp(str, "on")   == 0  ||
                               strcasecmp(str, "yes")  == 0  ||
                               strcasecmp(str, "true") == 0);
}


static EURLScheme s_ParseScheme(const char* str, size_t len)
{
    if (len == 5  &&  strncasecmp(str, "https", len) == 0)
        return eURL_Https;
    if (len == 4  &&  strncasecmp(str, "http",  len) == 0)
        return eURL_Http;
    if (len == 4  &&  strncasecmp(str, "file",  len) == 0)
        return eURL_File;
    if (len == 3  &&  strncasecmp(str, "ftp",   len) == 0)
        return eURL_Ftp;
    return eURL_Unspec;
}


static const char* s_Scheme(EURLScheme scheme)
{
    switch (scheme) {
    case eURL_Unspec:
        break;
    case eURL_Https:
        return "https";
    case eURL_Http:
        return "http";
    case eURL_File:
        return "file";
    case eURL_Ftp:
        return "ftp";
    }
    return 0;
}


/****************************************************************************
 * ConnNetInfo API
 */


extern SConnNetInfo* ConnNetInfo_Create(const char* service)
{
#define REG_VALUE(name, value, def_value)                               \
    ConnNetInfo_GetValue(service, name, value, sizeof(value), def_value)

    SConnNetInfo* info =
        (SConnNetInfo*) malloc(sizeof(*info)
                               + (service ? strlen(service) + 1 : 0));
    /* aux. storage */
    char   str[1024];
    size_t len;
    int    val;
    double dbl;
    char*  e;

    if (!info)
        return 0/*failure*/;

    /* client host: default */
    info->client_host[0] = '\0';

    /* scheme */
    REG_VALUE(REG_CONN_SCHEME, str, DEF_CONN_SCHEME);
    info->scheme = s_ParseScheme(str, strlen(str));

    /* username */
    REG_VALUE(REG_CONN_USER, info->user, DEF_CONN_USER);

    /* password */
    REG_VALUE(REG_CONN_PASS, info->pass, DEF_CONN_PASS);

    /* hostname */
    REG_VALUE(REG_CONN_HOST, info->host, DEF_CONN_HOST);

    /* port # */
    REG_VALUE(REG_CONN_PORT, str, DEF_CONN_PORT);
    errno = 0;
    if (*str  &&  (val = strtoul(str, &e, 10)) > 0  &&  !errno
        &&  !*e  &&  val < (1 << 16)) {
        info->port = val;
    } else
        info->port = 0/*default*/;

    /* path */
    REG_VALUE(REG_CONN_PATH, info->path, DEF_CONN_PATH);

    /* args */
    REG_VALUE(REG_CONN_ARGS, info->args, DEF_CONN_ARGS);

    /* request method */
    REG_VALUE(REG_CONN_REQ_METHOD, str, DEF_CONN_REQ_METHOD);
    if (!*str  ||  strcasecmp(str, "ANY") == 0)
        info->req_method = eReqMethod_Any;
    else if (strcasecmp(str, "POST") == 0)
        info->req_method = eReqMethod_Post;
    else if (strcasecmp(str, "GET") == 0)
        info->req_method = eReqMethod_Get;

    /* connection timeout */
    REG_VALUE(REG_CONN_TIMEOUT, str, 0);
    len = strlen(str);
    if (len > 2  &&  len < 9  &&  strncasecmp(str, "infinite", len) == 0) {
        info->timeout = 0;
    } else {
        info->timeout = &info->tmo;
        if (!*str  ||  (dbl = atof(str)) < 0.0)
            dbl = DEF_CONN_TIMEOUT;
        info->timeout->sec  = (unsigned int) dbl;
        info->timeout->usec = (unsigned int)
            ((dbl - info->timeout->sec) * 1000000.0);
    }

    /* max. # of attempts to establish connection */
    REG_VALUE(REG_CONN_MAX_TRY, str, 0);
    val = atoi(str);
    info->max_try = (unsigned short)(val > 0 ? val : DEF_CONN_MAX_TRY);

    /* HTTP proxy server? */
    REG_VALUE(REG_CONN_HTTP_PROXY_HOST, info->http_proxy_host,
              DEF_CONN_HTTP_PROXY_HOST);
    if (*info->http_proxy_host) {
        /* yes, use the specified HTTP proxy server */
        REG_VALUE(REG_CONN_HTTP_PROXY_PORT, str, DEF_CONN_HTTP_PROXY_PORT);
        errno = 0;
        if (*str  &&  (val = strtoul(str, &e, 10)) > 0
            &&  !errno  &&  !*e  &&  val < (1 << 16)) {
            info->http_proxy_port = val;
        } else
            info->http_proxy_port = 0/*default*/;
    } else
        info->http_proxy_port = 0;

    /* non-transparent CERN-like firewall proxy server? */
    REG_VALUE(REG_CONN_PROXY_HOST, info->proxy_host, DEF_CONN_PROXY_HOST);

    /* turn on debug printout? */
    REG_VALUE(REG_CONN_DEBUG_PRINTOUT, str, DEF_CONN_DEBUG_PRINTOUT);
    if (ConnNetInfo_Boolean(str)
        ||    (*str  &&   strcasecmp(str, "some") == 0)) {
        info->debug_printout = eDebugPrintout_Some;
    } else if (*str  &&  (strcasecmp(str, "all")  == 0  ||
                          strcasecmp(str, "data") == 0)) {
        info->debug_printout = eDebugPrintout_Data;
    } else
        info->debug_printout = eDebugPrintout_None;

    /* stateless client? */
    REG_VALUE(REG_CONN_STATELESS, str, DEF_CONN_STATELESS);
    info->stateless = ConnNetInfo_Boolean(str);

    /* firewall mode? */
    REG_VALUE(REG_CONN_FIREWALL, str, DEF_CONN_FIREWALL);
    info->firewall = ConnNetInfo_Boolean(str);

    /* prohibit the use of local load balancer? */
    REG_VALUE(REG_CONN_LB_DISABLE, str, DEF_CONN_LB_DISABLE);
    info->lb_disable = ConnNetInfo_Boolean(str);

    /* user header */
    REG_VALUE(REG_CONN_HTTP_USER_HEADER, str, DEF_CONN_HTTP_USER_HEADER);
    info->http_user_header = *str ? x_StrcatCRLF(NULL, str) : 0;

    /* default referer */
    ConnNetInfo_GetValue(0, REG_CONN_HTTP_REFERER, str, sizeof(str),
                         DEF_CONN_HTTP_REFERER);
    info->http_referer = *str ? strdup(str) : 0;

    /* not adjusted yet... */
    info->http_proxy_adjusted = 0/*false*/;

    /* store the service name for which this structure has been created */
    info->service = service ? strcpy((char*)info + sizeof(*info), service) : 0;

    /* done */
    return info;
#undef REG_VALUE
}


extern int/*bool*/ ConnNetInfo_AdjustForHttpProxy(SConnNetInfo* info)
{
    char   port[10];
    size_t hostlen;
    size_t portlen;
    size_t pathlen;
    size_t len;

    if (!info)
        return 0/*false*/;

    if (!*info->http_proxy_host  ||  info->http_proxy_adjusted)
        return 1/*true*/;

    if (info->scheme != eURL_Unspec  &&  info->scheme != eURL_Http) {
        CORE_LOGF_X(13, eLOG_Error,
                    ("[ConnNetInfo_AdjustForHttpProxy] "
                     " Cannot adjust for scheme \"%s\"",
                     s_Scheme(info->scheme)));
        assert(0);
        return 0/*false*/;
    }

    hostlen = strlen(info->host);
    portlen = info->port ? (size_t) sprintf(port, ":%hu", info->port) : 0;
    if (*info->path != '/')
        port[portlen++] = '/';
    pathlen = strlen(info->path);

    len = 7/*http://*/ + hostlen + portlen + pathlen + 1;

    if (len >= sizeof(info->path)) {
        CORE_LOGF_X(14, eLOG_Error,
                    ("[ConnNetInfo_AdjustForHttpProxy] "
                     "  Adjusted path too long (%lu)",
                     (unsigned long) len));
        assert(0);
        return 0/*false*/;
    }
    len -= pathlen + 1;
    memmove(info->path + len, info->path, pathlen + 1);
    len -= portlen;
    memcpy (info->path + len, port,       portlen);
    assert(len - hostlen == 7);
    memcpy (info->path + 7,   info->host, hostlen);
    memcpy (info->path,       "http://",  7);

    assert(sizeof(info->host) >= sizeof(info->http_proxy_host));
    strncpy0(info->host, info->http_proxy_host, sizeof(info->host) - 1);
    info->port = info->http_proxy_port;
    info->http_proxy_adjusted = 1/*true*/;
    info->scheme = eURL_Http;
    return 1/*true*/;
}


extern int/*bool*/ ConnNetInfo_ParseURL(SConnNetInfo* info, const char* url)
{
    const char *s, *a;
    size_t len;
    char* p;

    if (info->http_proxy_adjusted) {
        /* undo proxy adjustment */
        SConnNetInfo* temp = ConnNetInfo_Create(info->service);
        if (!ConnNetInfo_ParseURL(temp, info->path)) {
            ConnNetInfo_Destroy(temp);
            return 0/*failure*/;
        }
        memcpy(info->host, temp->host, sizeof(info->host));
        info->port = temp->port;
        memcpy(info->path, temp->path, sizeof(info->path));
        ConnNetInfo_Destroy(temp);
        info->http_proxy_adjusted = 0/*false*/;
    }

    /* "user:pass@host:port" first [any optional] */
    if ((s = strstr(url, "://")) != 0) {
        const char* h = s + 3; /*host starts here*/

        len = (size_t)(s - url);
        if ((info->scheme = s_ParseScheme(url, len)) == eURL_Unspec)
            return 0/*failure*/;

        /* find the end of host spec */
        len = strcspn(h, "/?#");
        s = h + len;

        /* username:password */
        if (!(a = (const char*) memchr(h, '@', len))) {
            info->user[0] = '\0';
            info->pass[0] = '\0';
        } else {
            size_t   ulen = (size_t)(a - h);
            const char* t = (const char*) memchr(h, ':', ulen);
            if (t)
                ulen = (size_t)(t - h);
            if (ulen < sizeof(info->user)) {
                memcpy(info->user, h, ulen);
                info->user[ulen] = '\0';
            } else {
                memcpy(info->user, h, sizeof(info->user) - 1);
                info->user[sizeof(info->user) - 1] = '\0';
            }
            if (t  &&  ulen) {
                /* take pass only if user non-empty */
                ulen = (size_t)(a - ++t);
                if (ulen < sizeof(info->pass)) {
                    memcpy(info->pass, t, ulen);
                    info->pass[ulen] = '\0';
                } else {
                    memcpy(info->pass, t, sizeof(info->pass) - 1);
                    info->pass[sizeof(info->pass) - 1] = '\0';
                }
            } else
                info->pass[0] = '\0';
            h = ++a;
            len = (size_t)(s - h);
        }

        /* host ends at "s" here */
        if (!(a = (const char*) memchr(h, ':', len))) {
            info->port = 0/*default*/;
            /*len remains unchanged*/
        } else if (isdigit(a[1])) {
            unsigned short port;
            int n;
            if (sscanf(a, ":%hu%n", &port, &n) < 1  ||  a + n != s  ||  !port)
                return 0/*failure*/;
            info->port = port;
            len = (size_t)(a - h);
        } else
            return 0/*failure*/;
        if (len < sizeof(info->host)) {
            memcpy(info->host, h, len);
            info->host[len] = '\0';
        } else {
            memcpy(info->host, h, sizeof(info->host) - 1);
            info->host[sizeof(info->host) - 1] = '\0';
        }
    } else
        s = url;

    a = s + strcspn(s, "?#");
    len = (size_t)(a - s);

    /* path (NB: can be relative); "len" holds the path length */
    if (s != url  ||  *s == '/'  ||  !(p = strrchr(info->path, '/'))) {
        /* absolute path */
        p = info->path;
        if (!len) {
            s = "/";   /* in case of an empty path we take the root '/' */
            len = 1;
        }
    } else
        p++;
    if (len < sizeof(info->path) - (size_t)(p - info->path)) {
        memcpy(p, s, len);
        p[len] = '\0';
    } else {
        memcpy(p, s, sizeof(info->path) - (size_t)(p - info->path) - 1);
        info->path[sizeof(info->path) - 1] = '\0';
    }

    /* arguments and fragment */
    if (*a) {
        if (*a == '#') {
            len = 0/*unused*/;
            s = a;
        } else if (!(s = strchr(++a/*NB: *a=='?'*/, '#'))) {
            len = strlen(a);
            s = a + len;
        } else
            len = (size_t)(s - a);
        /* "len" holds the length of new args("a", w/o the leading "?") */
        assert(!*s  ||  *s == '#');
        assert(*a != '?');

        if (*s) {
            /* if there is a new fragment, the entire args get overridden */
            if (!s[1]) {
                /* don't store the empty fragment # */
                len = (size_t)(s - a);
                if (len > sizeof(info->args) - 1)
                    len = sizeof(info->args) - 1;
            } else
                len = sizeof(info->args) - 1;
            strncpy0(info->args, a, len);
        } else if (!(p = strchr(info->args, '#'))
                   ||  len >= sizeof(info->args) - 1) {
            /* there is no new fragment and: either there is no old fragment,
             * or the old one will not fit -- in both cases, replace */
            strncpy0(info->args, a, sizeof(info->args) - 1);
        } else {
            /* there is no new fragment, but there was an old one -- keep it */
            size_t move = strlen(p);
            if (move > sizeof(info->args) - 1 - len)
                move = sizeof(info->args) - 1 - len;
            memmove(info->args + len, p, move);
            memcpy (info->args,       a, len);
            info->args[len + move] = '\0';
        }
    } else if ((p = strchr(info->args, '#')) != 0) {
        /* keep the old fragment, if any, but drop all args */
        memmove(info->args, p, strlen(p) + 1);
    }

    return 1/*success*/;
}


extern int/*bool*/ ConnNetInfo_SetUserHeader(SConnNetInfo* info,
                                             const char*   user_header)
{
    if (info->http_user_header)
        free((void*) info->http_user_header);
    if (!user_header  ||  !*user_header)
        info->http_user_header = 0;
    else if (!(info->http_user_header = x_StrcatCRLF(NULL, user_header)))
        return 0/*failure*/;
    return 1/*success*/;
}


extern int/*bool*/ ConnNetInfo_AppendUserHeader(SConnNetInfo* info,
                                                const char*   user_header)
{
    char* new_header;

    if (!info->http_user_header  ||  !*info->http_user_header)
        return ConnNetInfo_SetUserHeader(info, user_header);

    new_header = (char*) info->http_user_header;
    if (!(new_header = x_StrcatCRLF(new_header, user_header)))
        return 0/*failure*/;

    info->http_user_header = new_header;
    return 1/*success*/;
}


typedef enum {
    eUserHeaderOp_Delete,
    eUserHeaderOp_Extend,
    eUserHeaderOp_Override
} EUserHeaderOp;


static int/*bool*/ s_ModifyUserHeader(SConnNetInfo* info,
                                      const char*   user_header,
                                      EUserHeaderOp op)
{
    int/*bool*/ retval;
    char*  new_header;
    size_t newlinelen;
    size_t newhdrlen;
    char*  newline;
    size_t hdrlen;
    char*  hdr;

    if (!user_header || !(newhdrlen = strlen(user_header)))
        return 1/*success*/;

    if (!(hdr = (char*) info->http_user_header) || !(hdrlen = strlen(hdr))) {
        if (op == eUserHeaderOp_Delete)
            return 1/*success*/;
        if (!hdr && !(hdr = strdup("")))
            return 0/*failure*/;
        hdrlen = 0;
    }

    if (op != eUserHeaderOp_Delete) {
        if (!(new_header = (char*) malloc(newhdrlen + 1)))
            return 0/*failure*/;
        memcpy(new_header, user_header, newhdrlen + 1);
    } else
        new_header = (char*) user_header; /* we actually won't modify it! */

    retval = 1/*assume best: success*/;
    for (newline = new_header; *newline; newline += newlinelen) {
        char*  eol = strchr(newline, '\n');
        char*  eot = strchr(newline,  ':');
        int/*bool*/ used = 0;
        size_t newtaglen;
        char*  newtagval;
        size_t linelen;
        char*  line;
        size_t len;
        size_t l;

        newlinelen = (size_t)
            (eol ? eol - newline + 1 : new_header + newhdrlen - newline);
        if (!eot || eot >= newline + newlinelen ||
            !(newtaglen = (size_t)(eot - newline)))
            continue;

        newtagval = newline + newtaglen + 1;
        while (newtagval < newline + newlinelen) {
            if (isspace((unsigned char)(*newtagval)))
                newtagval++;
            else
                break;
        }
        switch (op) {
        case eUserHeaderOp_Delete:
            len = 0;
            break;
        case eUserHeaderOp_Extend:
            len = newlinelen - (size_t)(newtagval - newline);
            break;
        case eUserHeaderOp_Override:
            len = newtagval < newline + newlinelen ? newlinelen : 0;
            break;
        default:
            assert(0);
            retval = 0/*failure*/;
            len = 0;
            break;
        }

        for (line = hdr; *line; line += linelen) {
            size_t taglen;

            eol = strchr(line, '\n');
            eot = strchr(line,  ':');

            linelen = (size_t)(eol ? eol - line + 1 : hdr + hdrlen - line);
            if (!eot || eot >= line + linelen)
                continue;

            taglen = (size_t)(eot - line);
            if (newtaglen != taglen || strncasecmp(newline, line, taglen) != 0)
                continue;

            if (op == eUserHeaderOp_Extend) {
                if (linelen - taglen >= len
                    &&  strncasecmp(line + linelen - len, newtagval, len) == 0
                    &&  (linelen - taglen == len  ||
                         isspace((unsigned char) line[linelen - len - 1]))) {
                    len = 0;
                }
                l = linelen + len;
                if (len && linelen > 1 && line[linelen - 2] == '\r')
                    --l;
            } else
                l = len;
            if (l != linelen) {
                size_t off = (size_t)(line - hdr);
                if (l > linelen) {
                    char* temp = (char*)realloc(hdr, hdrlen + l - linelen + 1);
                    if (!temp) {
                        retval = 0/*failure*/;
                        continue;
                    }
                    hdr  = temp;
                    line = temp + off;
                }
                hdrlen -= linelen;
                memmove(line + l, line + linelen, hdrlen - off + 1);
                hdrlen += l;
            }

            if (len) {
                if (op == eUserHeaderOp_Extend) {
                    char* s = &line[l - len - 1];
                    *s++ = ' ';
                    memcpy(s, newtagval, len);
                } else
                    memcpy(line, newline, len);
                linelen = l;
                used = 1;
            } else if (op == eUserHeaderOp_Extend)
                used = 1;
        }

        if (op == eUserHeaderOp_Delete)
            continue;

        if (used || !len) {
            memmove(newline, newline + newlinelen,
                    newhdrlen - (size_t)(newline-new_header) - newlinelen + 1);
            newhdrlen -= newlinelen;
            newlinelen = 0;
        }
    }

    info->http_user_header = hdr;
    if (op != eUserHeaderOp_Delete) {
        if (!ConnNetInfo_AppendUserHeader(info, new_header))
            retval = 0/*failure*/;
        free(new_header);
    }
    return retval;
}


extern int/*bool*/ ConnNetInfo_OverrideUserHeader(SConnNetInfo* info,
                                                  const char*   header)
{
    return s_ModifyUserHeader(info, header, eUserHeaderOp_Override);
}


extern void ConnNetInfo_DeleteUserHeader(SConnNetInfo* info,
                                         const char*   header)
{
    verify(s_ModifyUserHeader(info, header, eUserHeaderOp_Delete));
}


extern int/*bool*/ ConnNetInfo_ExtendUserHeader(SConnNetInfo* info,
                                                const char*   header)
{
    return s_ModifyUserHeader(info, header, eUserHeaderOp_Extend);
}


extern int/*bool*/ ConnNetInfo_AppendArg(SConnNetInfo* info,
                                         const char*   arg,
                                         const char*   val)
{
    size_t len, used;

    if (!arg || !*arg)
        return 1/*success*/;

    used = strlen(info->args);
    len  = strlen(arg);
    
    if (used + (used ? 1/*&*/ : 0) + len +
        (val && *val ? 1/*=*/ + strlen(val) : 0) >= sizeof(info->args)) {
        return 0/*failure*/;
    }

    if (used)
        info->args[used++] = '&';
    strcpy(info->args + used, arg);
    if (val && *val) {
        used += len;
        info->args[used++] = '=';
        strcpy(info->args + used, val);
    }
    return 1/*success*/;
}


extern int/*bool*/ ConnNetInfo_PrependArg(SConnNetInfo* info,
                                          const char*   arg,
                                          const char*   val)
{
    size_t len, off, used;

    if (!arg || !*arg)
        return 1/*success*/;

    used = strlen(info->args);
    len  = strlen(arg);
    off  = len + (val && *val ? 1/*=*/ + strlen(val) : 0) + (used? 1/*&*/ : 0);

    if (off + used >= sizeof(info->args))
        return 0/*failure*/;

    if (used)
        memmove(info->args + off, info->args, used + 1);
    strcpy(info->args, arg);
    if (val && *val) {
        info->args[len++] = '=';
        strcpy(info->args + len, val);
    }
    if (used)
        info->args[off - 1] = '&';
    return 1/*success*/;
}


extern void ConnNetInfo_DeleteArg(SConnNetInfo* info,
                                  const char*   arg)
{
    size_t argnamelen;
    size_t arglen;
    char*  a;

    if (!arg || !(argnamelen = strcspn(arg, "=&")))
        return;
    for (a = info->args; *a; a += arglen) {
        if (*a == '&')
            a++;
        arglen = strcspn(a, "&");
        if (arglen < argnamelen  ||  strncasecmp(a, arg, argnamelen) != 0  ||
            (a[argnamelen] && a[argnamelen] != '=' && a[argnamelen] != '&'))
            continue;
        if (a[arglen]) {
            arglen++;     /* for intermediary args, eat '&' separator, too */
            memmove(a, a + arglen, strlen(a + arglen) + 1);
        } else if (a != info->args) {
            *--a = '\0';  /* last argument in a list: remove trailing '&' */
        } else {
            *a = '\0';    /* last and the only argument removed */
        }
        arglen = 0;
    }
}


extern void ConnNetInfo_DeleteAllArgs(SConnNetInfo* info,
                                      const char*   args)
{
    char* temp;
    char* arg;
    if (!args || !*args || !(temp = strdup(args))) {
        if (args && *args)
            *info->args = '\0';
        return;
    }
    arg = temp;
    while (*arg) {
        char* end = strchr(arg, '&');
        if (!end)
            end = arg + strlen(arg);
        else
            *end++ = '\0';
        ConnNetInfo_DeleteArg(info, arg);
        arg = end;
    }
    free(temp);
}


extern int/*bool*/ ConnNetInfo_PreOverrideArg(SConnNetInfo* info,
                                              const char*   arg,
                                              const char*   val)
{
    if (!arg || !*arg)
        return 1/*success*/;
    ConnNetInfo_DeleteAllArgs(info, arg);
    return ConnNetInfo_PrependArg(info, arg, val);
}


extern int/*bool*/ ConnNetInfo_PostOverrideArg(SConnNetInfo* info,
                                               const char*   arg,
                                               const char*   val)
{
    if (!arg || !*arg)
        return 1/*success*/;
    ConnNetInfo_DeleteAllArgs(info, arg);
    return ConnNetInfo_AppendArg(info, arg, val);
}


static int/*bool*/ s_IsSufficientAddress(const char* addr)
{
    const char*c;
    return (SOCK_isip(addr)  ||
            ((c = strchr(addr, '.'))  != 0  &&  c[1]  &&
             (c = strchr(c + 2, '.')) != 0  &&  c[1]));
}


static const char* s_ClientAddress(const char* client_host,
                                   int/*bool*/ local_host)
{
    const char* c = client_host;
    unsigned int ip;
    char addr[80];
    char* s;

    assert(client_host);
    strncpy0(addr, client_host, sizeof(addr) - 1);
    if (UTIL_NcbiLocalHostName(addr)  &&  (s = strdup(addr)) != 0)
        client_host = s;
    if (s_IsSufficientAddress(client_host)                          ||
        !(ip = *client_host  &&  !local_host
          ? SOCK_gethostbyname(client_host)
          : SOCK_GetLocalHostAddress(eDefault))                     ||
        SOCK_ntoa(ip, addr, sizeof(addr)) != 0                      ||
        !(s = (char*) malloc(strlen(client_host) + strlen(addr) + 3))) {
        return client_host;
    }
    sprintf(s, "%s(%s)", client_host, addr);
    if (c != client_host)
        free((void*) client_host);
    return s;
}


extern int/*bool*/ ConnNetInfo_SetupStandardArgs(SConnNetInfo* info,
                                                 const char*   service)
{
    static const char kService[]  = "service";
    static const char kAddress[]  = "address";
    static const char kPlatform[] = "platform";
    int/*bool*/ local_host;
    const char* arch;
    const char* addr;

    if (!info)
        return 0/*failed*/;
    /* Dispatcher CGI args (may sacrifice some if they don't fit altogether) */
    if (!(arch = CORE_GetPlatform())  ||  !*arch)
        ConnNetInfo_DeleteArg(info, kPlatform);
    else
        ConnNetInfo_PreOverrideArg(info, kPlatform, arch);
    local_host = !info->client_host[0];
    if (local_host  &&
        !SOCK_gethostbyaddr(0, info->client_host, sizeof(info->client_host))) {
        SOCK_gethostname(info->client_host, sizeof(info->client_host));
    }
    if (!(addr = s_ClientAddress(info->client_host, local_host))  ||  !*addr)
        ConnNetInfo_DeleteArg(info, kAddress);
    else
        ConnNetInfo_PreOverrideArg(info, kAddress, addr);
    if (addr != info->client_host)
        free((void*) addr);
    if (service  &&  *service) {
        if (!ConnNetInfo_PreOverrideArg(info, kService, service)) {
            ConnNetInfo_DeleteArg(info, kPlatform);
            if (!ConnNetInfo_PreOverrideArg(info, kService, service)) {
                ConnNetInfo_DeleteArg(info, kAddress);
                if (!ConnNetInfo_PreOverrideArg(info, kService, service))
                    return 0/*failed*/;
            }
        }
    }
    return 1/*succeeded*/;
}


extern SConnNetInfo* ConnNetInfo_Clone(const SConnNetInfo* info)
{
    SConnNetInfo* x_info;
    const char* service;

    if (!info)
        return 0;

    service = info->service;
    x_info = (SConnNetInfo*) malloc(sizeof(*x_info)
                                    + (service ? strlen(service) + 1 : 0));
    if (!x_info)
        return 0;

    memcpy(x_info, info, sizeof(*x_info));
    x_info->http_user_header = 0;
    x_info->http_referer = 0;

    if (info->timeout  &&  info->timeout != kDefaultTimeout) {
        x_info->tmo     = *info->timeout;
        x_info->timeout = &x_info->tmo;
    }
    if (info->http_user_header
        &&  !(x_info->http_user_header = strdup(info->http_user_header))) {
        ConnNetInfo_Destroy(x_info);
        return 0;
    }
    if (x_info->http_referer
        &&  !(x_info->http_referer = strdup(info->http_referer))) {
        ConnNetInfo_Destroy(x_info);
        return 0;
    }
    x_info->service
        = service ? strcpy((char*) x_info + sizeof(*x_info), service) : 0;
    return x_info;
}


static const char* s_SchemeStr(EURLScheme scheme, char buf[])
{
    switch (scheme) {
    case eURL_Ftp:
        return "ftp";
    case eURL_File:
        return "file";
    case eURL_Http:
        return "http";
    case eURL_Https:
        return "https";
    case eURL_Unspec:
        break;
    default:
        sprintf(buf, "?#0x%08X", (int) scheme);
        return buf;
    }
    return 0;
}

static const char* s_PortStr(unsigned short port, char buf[])
{
    if (port) {
        sprintf(buf, "%hu", port);
        return buf;
    }
    return "(default)";
}

static void s_SaveStringQuot(char* s, const char* name,
                             const char* str, int/*bool*/ quote)
{
    sprintf(s + strlen(s), "%-16.16s: %s%s%s\n", name,
            str && quote ? "\"" : "",
            str          ? str  : "NULL",
            str && quote ? "\"" : "");
}

static void s_SaveString(char* s, const char* name, const char* str)
{
    s_SaveStringQuot(s, name, str, 1);
}

static void s_SaveKeyval(char* s, const char* name, const char* str)
{
    s_SaveStringQuot(s, name, str, 0);
}

static void s_SaveULong(char* s, const char* name, unsigned long lll)
{
    sprintf(s + strlen(s), "%-16.16s: %lu\n", name, lll);
}

static void s_SaveBool(char* s, const char* name, int/*bool*/ bbb)
{
    sprintf(s + strlen(s), "%-16.16s: %s\n", name, bbb ? "TRUE" : "FALSE");
}

static void s_SaveUserHeader(char* s, const char* name,
                             const char* uh, size_t uhlen)
{
    s += strlen(s);
    s += sprintf(s, "%-16.16s: ", name);
    if (uh) {
        *s++ = '"';
        memcpy(UTIL_PrintableString(uh, uhlen, s, 0/*reduce*/), "\"\n", 3);
    } else
        memcpy(s, "NULL\n", 6);
}

extern void ConnNetInfo_LogEx(const SConnNetInfo* info, ELOG_Level sev, LOG lg)
{
    char   scheme[32];
    char   port[16];
    size_t uhlen;
    size_t len;
    char*  s;

    if (!lg) {
        if (sev == eLOG_Fatal)
            abort();
        return;
    }

    if (!info) {
        LOG_Write(lg, NCBI_C_ERRCODE_X, 10, sev, 0, 0, 0,
                  "ConnNetInfo_Log: NULL info", 0, 0);
        return;
    }

    uhlen = info->http_user_header ? strlen(info->http_user_header) : 0;
       
    len = sizeof(*info) + 1024/*slack for all labels*/
        + UTIL_PrintableStringSize(info->http_user_header, uhlen)
        + (info->http_referer ? strlen(info->http_referer) : 0)
        + (info->service ? strlen(info->service) : 0);

    if (!(s = (char*) malloc(len))) {
        LOG_WRITE(lg, NCBI_C_ERRCODE_X, 11,
                  sev == eLOG_Fatal ? eLOG_Fatal : eLOG_Error,
                  "ConnNetInfo_Log: Cannot allocate temporary buffer");
        return;
    }

    strcpy(s, "ConnNetInfo_Log\n"
           "#################### [BEGIN] SConnNetInfo:\n");
    s_SaveString    (s, "service",         info->service);
    if (*info->client_host)
        s_SaveString(s, "client_host",     info->client_host);
    else
        s_SaveKeyval(s, "client_host",     "(default)");
    s_SaveString    (s, "scheme",          s_SchemeStr(info->scheme, scheme));
    s_SaveString    (s, "user",            info->user);
    if (*info->pass)
        s_SaveKeyval(s, "pass",           *info->user ? "(set)" : "(ignored)");
    else
        s_SaveString(s, "pass",            info->pass);
    s_SaveString    (s, "host",            info->host);
    s_SaveKeyval    (s, "port",            s_PortStr(info->port, port));
    s_SaveString    (s, "path",            info->path);
    s_SaveString    (s, "args",            info->args);
    s_SaveKeyval    (s, "req_method",     (info->req_method == eReqMethod_Any
                                           ? "ANY"
                                           : (info->req_method
                                              == eReqMethod_Get
                                              ? "GET"
                                              : (info->req_method
                                                 == eReqMethod_Post
                                                 ? "POST" : "(unknown)"))));
    if (info->timeout) {
        s_SaveULong (s, "timeout(sec)",    info->timeout->sec);
        s_SaveULong (s, "timeout(usec)",   info->timeout->usec);
    } else
        s_SaveKeyval(s, "timeout",         "INFINITE");
    s_SaveULong     (s, "max_try",         info->max_try);
    s_SaveString    (s, "http_proxy_host", info->http_proxy_host);
    s_SaveKeyval    (s, "http_proxy_port", s_PortStr(info->http_proxy_port,
                                                     port));
    s_SaveString    (s, "proxy_host",      info->proxy_host);
    s_SaveKeyval    (s, "debug_printout", (info->debug_printout
                                           == eDebugPrintout_None
                                           ? "NONE"
                                           : (info->debug_printout
                                              == eDebugPrintout_Some
                                              ? "SOME"
                                              : (info->debug_printout
                                                 == eDebugPrintout_Data
                                                 ? "DATA" : "(unknown)"))));
    s_SaveBool      (s, "stateless",       info->stateless);
    s_SaveBool      (s, "firewall",        info->firewall);
    s_SaveBool      (s, "lb_disable",      info->lb_disable);
    s_SaveUserHeader(s, "http_user_header",info->http_user_header, uhlen);
    s_SaveString    (s, "http_referer",    info->http_referer);
    s_SaveBool      (s, "[proxy_adjusted]",info->http_proxy_adjusted);
    strcat(s, "#################### [END] SConnNetInfo\n");

    assert(strlen(s) < len);
    LOG_Write(lg, NCBI_C_ERRCODE_X, 12, sev, 0, 0, 0, s, 0, 0);
    free(s);
}


/*FIXME:  To remove*/
#undef ConnNetInfo_Log
extern void ConnNetInfo_Log(const SConnNetInfo* info, LOG lg)
{
    ConnNetInfo_LogEx(info, eLOG_Trace, lg);
}


extern char* ConnNetInfo_URL(const SConnNetInfo* info)
{
    const char* scheme;
    size_t      len;
    char*       url;

    if (!info  ||  info->scheme == eURL_Unspec
        ||  !(scheme = s_Scheme(info->scheme))) {
        return 0;
    }

    len = strlen(scheme) + 3/*"://"*/ + strlen(info->host)
        + (*info->user ? strlen(info->user) + strlen(info->pass) + 2 : 0)
        + (info->port ? 6/*:port*/ : 0)
        + (info->http_proxy_adjusted ? 2 : 0)
        + strlen(info->path)
        + (*info->args ? strlen(info->args) + 3 : 2);
    url = (char*) malloc(len);

    if (url) {
        len = (size_t) sprintf(url, "%s://%s", scheme,
                               *info->user ? info->user : info->host);
        if (*info->user)
            len += sprintf(url + len, ":%s@%s", info->pass, info->host);
        if (info->port)
            len += sprintf(url + len, ":%hu", info->port);
        sprintf(url + len, "%s%s%s%s%s", info->http_proxy_adjusted
                ? "(" : *info->path != '/' ? "/" : "", info->path,
                &"?"[!*info->args  ||  *info->args == '#'], info->args,
                ")" + !info->http_proxy_adjusted);
    }
    return url;
}


extern void ConnNetInfo_Destroy(SConnNetInfo* info)
{
    if (!info)
        return;
    ConnNetInfo_SetUserHeader(info, 0);
    if (info->http_referer) {
        free((void*) info->http_referer);
        info->http_referer = 0;
    }
    free(info);
}



/****************************************************************************
 * URL_Connect
 */


extern EIO_Status URL_ConnectEx
(const char*     host,
 unsigned short  port,
 const char*     path,
 const char*     args,
 EReqMethod      req_method,
 size_t          content_length,
 const STimeout* o_timeout,
 const STimeout* rw_timeout,
 const char*     user_hdr,
 int/*bool*/     encode_args,
 TSOCK_Flags     flags,
 SOCK*           sock)
{
    static const char X_REQ_Q[] = "?";
    static const char X_REQ_E[] = " HTTP/1.0\r\n";
    static const char X_HOST[]  = "Host: ";

    BUF         buf;
    char*       hdr;
    const char* temp;
    EIO_Status  status;
    size_t      hdr_len;
    char        hdr_buf[80];
    int         add_hdr = 1;
    size_t      args_len = args  &&  *args ? strcspn(args, "#") : 0;
    size_t      user_hdr_len = user_hdr  &&  *user_hdr ? strlen(user_hdr) : 0;
    const char* x_req_method; /* "POST "/"GET " */

    /* check the args */
    if (!sock  ||  !host  ||  !*host  ||
        (user_hdr  &&  *user_hdr  &&  user_hdr[user_hdr_len - 1] != '\n')) {
        CORE_LOG_X(2, eLOG_Error, "[URL_Connect]  Bad arguments");
        assert(0);
        return eIO_InvalidArg;
    }

    /* select request method and its verbal representation */
    if (req_method == eReqMethod_Any) {
        req_method = content_length ? eReqMethod_Post : eReqMethod_Get;
    } else if (req_method == eReqMethod_Get  &&  content_length) {
        CORE_LOG_X(3, eLOG_Warning,
                   "[URL_Connect]  Content length ignored with method GET");
        content_length = 0;
    }
    switch (req_method) {
    case eReqMethod_Post:
        x_req_method = "POST ";
        break;
    case eReqMethod_Get:
        x_req_method = "GET ";
        break;
    default:
        CORE_LOGF_X(4, eLOG_Error,
                    ("[URL_Connect]  Unrecognized request method (#%u)",
                     (unsigned int) req_method));
        assert(0);
        return eIO_InvalidArg;
    }

    for (temp = user_hdr;  temp  &&  *temp;  temp = strchr(temp, '\n')) {
        if (temp != user_hdr)
            temp++;
        if (strncasecmp(temp, X_HOST, sizeof(X_HOST) - 2) == 0) {
            add_hdr = 0;
            break;
        }
    }

    /* URL-encode "args", if any specified */
    if (args_len) {
        if (encode_args) {
            size_t rd_len, wr_len;
            size_t size = 3 * args_len;
            char* x_args = (char*) malloc(size);
            if (!x_args) {
                CORE_LOGF_ERRNO_X(8, eLOG_Error, errno,
                                  ("[URL_Connect]  Out of memory (%lu)",
                                   (unsigned long) size));
                return eIO_Unknown;
            }
            URL_Encode(args, args_len, &rd_len, x_args, size, &wr_len);
            assert(args_len == rd_len);
            args_len = wr_len;
            temp = x_args;
        } else
            temp = args;
    } else
        temp = 0;

    if (!port) {
        hdr_len = 0;
        port = flags & fSOCK_Secure ? 443 : 80;
    } else
        hdr_len = (size_t)(add_hdr ? sprintf(hdr_buf, ":%hu", port) : 0);

    buf = 0;
    errno = 0;
    /* compose HTTP header */
    if (/* {POST|GET} <path>?<args> HTTP/1.0\r\n */
        !BUF_Write(&buf, x_req_method, strlen(x_req_method))  ||
        !BUF_Write(&buf, path,         strlen(path))          ||
        (args_len
         &&  (!BUF_Write(&buf, X_REQ_Q, sizeof(X_REQ_Q) - 1)  ||
              !BUF_Write(&buf, temp,    args_len)))           ||
        !BUF_Write      (&buf, X_REQ_E, sizeof(X_REQ_E) - 1)  ||

        (add_hdr
         /* Host: host[:port]\r\n */
         &&  (!BUF_Write(&buf, X_HOST,  sizeof(X_HOST) - 1)   ||
              !BUF_Write(&buf, host,    strlen(host))         ||
              !BUF_Write(&buf, hdr_buf, hdr_len)              ||
              !BUF_Write(&buf, "\r\n",  2)))                  ||

        /* Content-Length: <content_length>\r\n */
        (req_method != eReqMethod_Get
         &&  ((add_hdr =
               sprintf(hdr_buf, "Content-Length: %lu\r\n",
                       (unsigned long) content_length)) < 0   ||
              !BUF_Write(&buf, hdr_buf, (size_t) add_hdr)))    ||

        /* <user_header> */
        (user_hdr_len
         &&  !BUF_Write(&buf, user_hdr, user_hdr_len))        ||

        /* header separator */
        !BUF_Write(&buf, "\r\n", 2)) {
        int x_errno = errno;
        CORE_LOGF_ERRNO_X(5, eLOG_Error, x_errno,
                          ("[URL_Connect]  Error building HTTP header for"
                           " %s:%hu", host, port));
        BUF_Destroy(buf);
        if (temp  &&  temp != args)
            free((void*) temp);
        return eIO_Unknown;
    }
    if (temp  &&  temp != args)
        free((void*) temp);

    if (!(hdr = (char*) malloc(hdr_len = BUF_Size(buf)))
        ||  BUF_Read(buf, hdr, hdr_len) != hdr_len) {
        int x_errno = errno;
        CORE_LOGF_ERRNO_X(6, eLOG_Error, x_errno,
                          ("[URL_Connect]  Error storing HTTP header for"
                           " %s:%hu", host, port));
        if (hdr)
            free(hdr);
        BUF_Destroy(buf);
        return eIO_Unknown;
    }
    BUF_Destroy(buf);

    /* connect to HTTPD */
    status = SOCK_CreateEx(host, port, o_timeout, sock, hdr, hdr_len, flags);
    free(hdr);

    if (status != eIO_Success) {
        char temp[80];
        assert(!*sock);
        if (status == eIO_Timeout  &&  o_timeout) {
            sprintf(temp, "[%u.%06u]",
                    (unsigned int)(o_timeout->sec + o_timeout->usec/1000000),
                    (unsigned int)                 (o_timeout->usec%1000000));
        } else
            *temp = '\0';
        CORE_LOGF_X(7, eLOG_Error,
                    ("[URL_Connect]  Socket connect to %s:%hu failed: %s%s",
                     host, port, IO_StatusStr(status), temp));
    } else
        verify(SOCK_SetTimeout(*sock, eIO_ReadWrite, rw_timeout)==eIO_Success);
    return status;
}


extern SOCK URL_Connect
(const char*     host,
 unsigned short  port,
 const char*     path,
 const char*     args,
 EReqMethod      req_method,
 size_t          content_length,
 const STimeout* c_timeout,
 const STimeout* rw_timeout,
 const char*     user_hdr,
 int/*bool*/     encode_args,
 TSOCK_Flags     flags)
{
    SOCK sock;
    EIO_Status st = URL_ConnectEx(host, port, path, args,
                                  req_method, content_length,
                                  c_timeout, rw_timeout,
                                  user_hdr, encode_args, flags, &sock);
    return st == eIO_Success ? sock : 0;
}



/****************************************************************************
 * StripToPattern()
 */


typedef EIO_Status (*FDoIO)
     (void*     stream,
      void*     buf,
      size_t    size,
      size_t*   n_read,
      EIO_Event what     /* eIO_Read | eIO_Write (to pushback) */
      );

static EIO_Status s_StripToPattern
(void*       stream,
 FDoIO       io_func,
 const void* pattern,
 size_t      pattern_size,
 BUF*        buf,
 size_t*     n_discarded)
{
    EIO_Status status;
    char*      buffer;
    size_t     buffer_size;
    size_t     n_read = 0;

    /* check args */
    if ( n_discarded )
        *n_discarded = 0;
    if (!stream  ||  (pattern != 0) != (pattern_size != 0))
        return eIO_InvalidArg;

    /* allocate a temporary read buffer */
    buffer_size = 2 * pattern_size;
    if (buffer_size < 4096)
        buffer_size = 4096;
    if ( !(buffer = (char*) malloc(buffer_size)) )
        return eIO_Unknown;

    if ( !pattern ) {
        /* read/discard until EOF */
        do {
            status = io_func(stream, buffer, buffer_size, &n_read, eIO_Read);
            if ( buf )
                BUF_Write(buf, buffer, n_read);
            if ( n_discarded )
                *n_discarded += n_read;
        } while (status == eIO_Success);
    } else {
        for (;;) {
            /* read; search for the pattern; store the discarded data */
            size_t x_read, n_stored;

            assert(n_read < pattern_size);
            status = io_func(stream, buffer + n_read, buffer_size - n_read,
                             &x_read, eIO_Read);
            if ( !x_read ) {
                assert(status != eIO_Success);
                break; /*error*/
            }
            n_stored = n_read + x_read;

            if (n_stored >= pattern_size) {
                /* search for the pattern */
                size_t n_check = n_stored - pattern_size + 1;
                const char* b;
                for (b = buffer;  n_check;  b++, n_check--) {
                    if (*b != *((const char*) pattern))
                        continue;
                    if (memcmp(b, pattern, pattern_size) == 0)
                        break; /*found*/
                }
                /* pattern found */
                if ( n_check ) {
                    size_t x_discarded = (size_t)(b - buffer) + pattern_size;
                    if ( buf )
                        BUF_Write(buf, buffer + n_read, x_discarded - n_read);
                    if ( n_discarded )
                        *n_discarded += x_discarded;
                    /* return unused portion to the stream */
                    status = io_func(stream, buffer + x_discarded,
                                     n_stored - x_discarded, 0, eIO_Write);
                    break; /*finished*/
                }
            }

            /* pattern not found yet */
            if ( buf )
                BUF_Write(buf, buffer + n_read, x_read);
            if ( n_discarded )
                *n_discarded += x_read;
            n_read = n_stored;

            if (n_read > pattern_size) {
                size_t n_cut = n_read - pattern_size + 1;
                n_read = pattern_size - 1;
                memmove(buffer, buffer + n_cut, n_read);
            }
        }
    }

    /* cleanup & exit */
    free(buffer);
    return status;
}


static EIO_Status s_CONN_IO
(void*     stream,
 void*     buf,
 size_t    size,
 size_t*   n_read,
 EIO_Event what)
{
    switch (what) {
    case eIO_Read:
        return CONN_Read((CONN) stream, buf, size, n_read, eIO_ReadPlain);
    case eIO_Write:
        assert(stream);
        return CONN_PushBack((CONN) stream, buf, size);
    default:
        break;
    }
    return eIO_InvalidArg;
}

extern EIO_Status CONN_StripToPattern
(CONN        conn,
 const void* pattern,
 size_t      pattern_size,
 BUF*        buf,
 size_t*     n_discarded)
{
    return s_StripToPattern
        (conn, s_CONN_IO, pattern, pattern_size, buf, n_discarded);
}


static EIO_Status s_SOCK_IO
(void*     stream,
 void*     buf,
 size_t    size,
 size_t*   n_read,
 EIO_Event what)
{
    switch (what) {
    case eIO_Read:
        return SOCK_Read((SOCK) stream, buf, size, n_read, eIO_ReadPlain);
    case eIO_Write:
        return SOCK_PushBack((SOCK) stream, buf, size);
    default:
        break;
    }
    return eIO_InvalidArg;
}

extern EIO_Status SOCK_StripToPattern
(SOCK        sock,
 const void* pattern,
 size_t      pattern_size,
 BUF*        buf,
 size_t*     n_discarded)
{
    return s_StripToPattern
        (sock, s_SOCK_IO, pattern, pattern_size, buf, n_discarded);
}


static EIO_Status s_BUF_IO
(void*     stream,
 void*     buf,
 size_t    size,
 size_t*   n_read,
 EIO_Event what)
{
    BUF b;
    switch (what) {
    case eIO_Read:
        *n_read = BUF_Read((BUF) stream, buf, size);
        return *n_read ? eIO_Success : eIO_Closed;
    case eIO_Write:
        assert(stream);
        b = (BUF) stream;
        return BUF_PushBack(&b, buf, size) ? eIO_Success : eIO_Unknown;
    default:
        break;
    }
    return eIO_InvalidArg;
}

extern EIO_Status BUF_StripToPattern
(BUF         buffer,
 const void* pattern,
 size_t      pattern_size,
 BUF*        buf,
 size_t*     n_discarded)
{
    return s_StripToPattern
        (buffer, s_BUF_IO, pattern, pattern_size, buf, n_discarded);
}




/****************************************************************************
 * URL- Encoding/Decoding
 */


/* Return integer (0..15) corresponding to the "ch" as a hex digit
 * Return -1 on error
 */
static int s_HexChar(char ch)
{
    unsigned int rc = ch - '0';
    if (rc <= 9)
        return rc;
    rc = (ch | ' ') - 'a';
    return rc <= 5 ? (int) rc + 10 : -1;
}


/* The URL-encoding table
 */
static const char s_EncodeTable[256][4] = {
    "%00", "%01", "%02", "%03", "%04", "%05", "%06", "%07",
    "%08", "%09", "%0A", "%0B", "%0C", "%0D", "%0E", "%0F",
    "%10", "%11", "%12", "%13", "%14", "%15", "%16", "%17",
    "%18", "%19", "%1A", "%1B", "%1C", "%1D", "%1E", "%1F",
    "+",   "!",   "%22", "%23", "$",   "%25", "%26", "'",
    "(",   ")",   "*",   "%2B", ",",   "-",   ".",   "%2F",
    "0",   "1",   "2",   "3",   "4",   "5",   "6",   "7",
    "8",   "9",   "%3A", "%3B", "%3C", "%3D", "%3E", "%3F",
    "%40", "A",   "B",   "C",   "D",   "E",   "F",   "G",
    "H",   "I",   "J",   "K",   "L",   "M",   "N",   "O",
    "P",   "Q",   "R",   "S",   "T",   "U",   "V",   "W",
    "X",   "Y",   "Z",   "%5B", "%5C", "%5D", "%5E", "_",
    "%60", "a",   "b",   "c",   "d",   "e",   "f",   "g",
    "h",   "i",   "j",   "k",   "l",   "m",   "n",   "o",
    "p",   "q",   "r",   "s",   "t",   "u",   "v",   "w",
    "x",   "y",   "z",   "%7B", "%7C", "%7D", "%7E", "%7F",
    "%80", "%81", "%82", "%83", "%84", "%85", "%86", "%87",
    "%88", "%89", "%8A", "%8B", "%8C", "%8D", "%8E", "%8F",
    "%90", "%91", "%92", "%93", "%94", "%95", "%96", "%97",
    "%98", "%99", "%9A", "%9B", "%9C", "%9D", "%9E", "%9F",
    "%A0", "%A1", "%A2", "%A3", "%A4", "%A5", "%A6", "%A7",
    "%A8", "%A9", "%AA", "%AB", "%AC", "%AD", "%AE", "%AF",
    "%B0", "%B1", "%B2", "%B3", "%B4", "%B5", "%B6", "%B7",
    "%B8", "%B9", "%BA", "%BB", "%BC", "%BD", "%BE", "%BF",
    "%C0", "%C1", "%C2", "%C3", "%C4", "%C5", "%C6", "%C7",
    "%C8", "%C9", "%CA", "%CB", "%CC", "%CD", "%CE", "%CF",
    "%D0", "%D1", "%D2", "%D3", "%D4", "%D5", "%D6", "%D7",
    "%D8", "%D9", "%DA", "%DB", "%DC", "%DD", "%DE", "%DF",
    "%E0", "%E1", "%E2", "%E3", "%E4", "%E5", "%E6", "%E7",
    "%E8", "%E9", "%EA", "%EB", "%EC", "%ED", "%EE", "%EF",
    "%F0", "%F1", "%F2", "%F3", "%F4", "%F5", "%F6", "%F7",
    "%F8", "%F9", "%FA", "%FB", "%FC", "%FD", "%FE", "%FF"
};

#define VALID_URL_SYMBOL(ch)  (s_EncodeTable[(unsigned char)ch][0] != '%')


extern int/*bool*/ URL_DecodeEx
(const void* src_buf,
 size_t      src_size,
 size_t*     src_read,
 void*       dst_buf,
 size_t      dst_size,
 size_t*     dst_written,
 const char* allow_symbols)
{
    unsigned char* src = (unsigned char*) src_buf;
    unsigned char* dst = (unsigned char*) dst_buf;

    *src_read    = 0;
    *dst_written = 0;
    if (!src_size  ||  !dst_size)
        return 1/*true*/;
    if (!src  ||  !dst)
        return 0/*false*/;

    for ( ;  *src_read != src_size  &&  *dst_written != dst_size;
          (*src_read)++, (*dst_written)++, src++, dst++) {
        switch ( *src ) {
        case '%': {
            int i1, i2;
            if (*src_read + 2 > src_size)
                return 1/*true*/;
            if ((i1 = s_HexChar(*(++src))) == -1)
                return *dst_written ? 1/*true*/ : 0/*false*/;
            if (*src_read + 3 > src_size)
                return 1/*true*/;
            if ((i2 = s_HexChar(*(++src))) == -1)
                return *dst_written ? 1/*true*/ : 0/*false*/;

            *dst = (unsigned char)((i1 << 4) + i2);
            *src_read += 2;
            break;
        }
        case '+': {
            *dst = ' ';
            break;
        }
        default: {
            if (VALID_URL_SYMBOL(*src)  ||
                (allow_symbols  &&  strchr(allow_symbols, *src)))
                *dst = *src;
            else
                return *dst_written ? 1/*true*/ : 0/*false*/;
        }
        }/*switch*/
    }

    assert(src == (unsigned char*) src_buf + *src_read   );
    assert(dst == (unsigned char*) dst_buf + *dst_written);
    return 1/*true*/;
}


extern int/*bool*/ URL_Decode
(const void* src_buf,
 size_t      src_size,
 size_t*     src_read,
 void*       dst_buf,
 size_t      dst_size,
 size_t*     dst_written)
{
    return URL_DecodeEx
        (src_buf, src_size, src_read, dst_buf, dst_size, dst_written, 0);
}


extern void URL_Encode
(const void* src_buf,
 size_t      src_size,
 size_t*     src_read,
 void*       dst_buf,
 size_t      dst_size,
 size_t*     dst_written)
{
    unsigned char* src = (unsigned char*) src_buf;
    unsigned char* dst = (unsigned char*) dst_buf;

    *src_read    = 0;
    *dst_written = 0;
    if (!src_size  ||  !dst_size  ||  !dst  ||  !src)
        return;

    for ( ;  *src_read != src_size  &&  *dst_written != dst_size;
          (*src_read)++, (*dst_written)++, src++, dst++) {
        const char* subst = s_EncodeTable[*src];
        if (*subst != '%') {
            *dst = *subst;
        } else if (*dst_written < dst_size - 2) {
            *dst = '%';
            *(++dst) = *(++subst);
            *(++dst) = *(++subst);
            *dst_written += 2;
        } else {
            return;
        }
    }
    assert(src == (unsigned char*) src_buf + *src_read   );
    assert(dst == (unsigned char*) dst_buf + *dst_written);
}


/****************************************************************************
 * NCBI-specific MIME content type and sub-types
 */


static const char* s_MIME_Type[eMIME_T_Unknown+1] = {
    "x-ncbi-data",
    "text",
    "application",
    "unknown"
};

static const char* s_MIME_SubType[eMIME_Unknown+1] = {
    "x-dispatch",
    "x-asn-text",
    "x-asn-binary",
    "x-fasta",
    "x-www-form",
    "html",
    "plain",
    "xml",
    "xml+soap",
    "octet-stream",
    "x-unknown"
};

static const char* s_MIME_Encoding[eENCOD_Unknown+1] = {
    "",
    "urlencoded",
    "encoded"
};


extern char* MIME_ComposeContentTypeEx
(EMIME_Type     type,
 EMIME_SubType  subtype,
 EMIME_Encoding encoding,
 char*          buf,
 size_t         buflen)
{
    static const char s_ContentType[] = "Content-Type: ";
    const char* x_type;
    const char* x_subtype;
    const char* x_encoding;
    char        x_buf[MAX_CONTENT_TYPE_LEN];

    assert(buf  &&  buflen);

    if (type == eMIME_T_Undefined  ||  subtype == eMIME_Undefined)
        return 0;
    if (type >= eMIME_T_Unknown)
        type  = eMIME_T_Unknown;
    if (subtype >= eMIME_Unknown)
        subtype  = eMIME_Unknown;
    if (encoding >= eENCOD_Unknown)
        encoding  = eENCOD_Unknown;

    x_type     = s_MIME_Type    [type];
    x_subtype  = s_MIME_SubType [subtype];
    x_encoding = s_MIME_Encoding[encoding];

    if ( *x_encoding ) {
        assert(sizeof(s_ContentType) + strlen(x_type) + strlen(x_subtype)
               + strlen(x_encoding) + 4 < MAX_CONTENT_TYPE_LEN);
        sprintf(x_buf, "%s%s/%s-%s\r\n",
                s_ContentType, x_type, x_subtype, x_encoding);
    } else {
        assert(sizeof(s_ContentType) + strlen(x_type) + strlen(x_subtype)
               + 3 < MAX_CONTENT_TYPE_LEN);
        sprintf(x_buf, "%s%s/%s\r\n", s_ContentType, x_type, x_subtype);
    }
    assert(strlen(x_buf) < sizeof(x_buf));
    assert(strlen(x_buf) < buflen);
    strncpy0(buf, x_buf, buflen - 1);
    return buf;
}


extern int/*bool*/ MIME_ParseContentTypeEx
(const char*     str,
 EMIME_Type*     type,
 EMIME_SubType*  subtype,
 EMIME_Encoding* encoding)
{
    char*  x_buf;
    size_t x_size;
    char*  x_type;
    char*  x_subtype;
    int    i;

    if ( type )
        *type = eMIME_T_Undefined;
    if ( subtype )
        *subtype = eMIME_Undefined;
    if ( encoding )
        *encoding = eENCOD_None;

    x_size = str  &&  *str ? strlen(str) + 1 : 0;
    if (!x_size)
        return 0/*false*/;

    if (!(x_buf = (char*) malloc(x_size << 1)))
        return 0/*false*/;
    x_type = x_buf + x_size;

    strlwr(strcpy(x_buf, str));

    if ((sscanf(x_buf, " content-type: %s ", x_type) != 1  &&
         sscanf(x_buf, " %s ", x_type) != 1)  ||
        (x_subtype = strchr(x_type, '/')) == 0) {
        free(x_buf);
        return 0/*false*/;
    }
    *x_subtype++ = '\0';
    x_size = strlen(x_subtype);

    if ( type ) {
        for (i = 0;  i < (int) eMIME_T_Unknown;  i++) {
            if (strcmp(x_type, s_MIME_Type[i]) == 0)
                break;
        }
        *type = (EMIME_Type) i;
    }

    for (i = 1;  i <= (int) eENCOD_Unknown;  i++) {
        size_t len = strlen(s_MIME_Encoding[i]);
        if (len < x_size) {
            char* x_encoding = x_subtype + x_size - len;
            if (x_encoding[-1] == '-'
                &&  strcmp(x_encoding, s_MIME_Encoding[i]) == 0) {
                if ( encoding ) {
                    *encoding = (i == (int) eENCOD_Unknown
                                 ? eENCOD_None : (EMIME_Encoding) i);
                }
                x_encoding[-1] = '\0';
                break;
            }
        }
    }

    if ( subtype ) {
        for (i = 0;  i < (int) eMIME_Unknown;  i++) {
            if (strcmp(x_subtype, s_MIME_SubType[i]) == 0)
                break;
        }
        *subtype = (EMIME_SubType) i;
    }

    free(x_buf);
    return 1/*true*/;
}


/* DEPRECATED and scheduled for removal */
extern char* MIME_ComposeContentType
(EMIME_SubType  subtype,
 EMIME_Encoding encoding,
 char*          buf,
 size_t         buflen)
{
    return MIME_ComposeContentTypeEx(eMIME_T_NcbiData,
                                     subtype, encoding, buf, buflen);
}


/* DEPRECATED and scheduled for removal */
extern int/*bool*/ MIME_ParseContentType
(const char*     str,
 EMIME_SubType*  subtype,
 EMIME_Encoding* encoding)
{
    EMIME_Type type;
    if ( !MIME_ParseContentTypeEx(str, &type, subtype, encoding) )
        return 0/*false*/;

    if (type != eMIME_T_NcbiData) {
        if ( subtype )
            *subtype  = eMIME_Unknown;
        if ( encoding )
            *encoding = eENCOD_Unknown;
        return 0/*false*/;
    }

    return 1/*true*/;
}
