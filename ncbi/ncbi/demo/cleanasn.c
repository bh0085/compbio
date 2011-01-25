/*   cleanasn.c
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*            National Center for Biotechnology Information (NCBI)
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government do not place any restriction on its use or reproduction.
*  We would, however, appreciate having the NCBI and the author cited in
*  any work or product based on this material
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
* ===========================================================================
*
* File Name:  cleanasn.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   10/19/99
*
* $Revision: 6.143 $
*
* File Description:
*
* Modifications:
* --------------------------------------------------------------------------
* Date     Name        Description of modification
* -------  ----------  -----------------------------------------------------
*
*
* ==========================================================================
*/

#include <ncbi.h>
#include <objall.h>
#include <objsset.h>
#include <objfdef.h>
#include <objsub.h>
#include <sequtil.h>
#include <gather.h>
#include <sqnutils.h>
#include <explore.h>
#include <tofasta.h>
#include <toasn3.h>
#include <toporg.h>
#include <subutil.h>
#include <asn2gnbk.h>
#include <pmfapi.h>
#include <tax3api.h>
#include <asn2gnbi.h>
#include <ent2api.h>
#include <gbftdef.h>
#ifdef INTERNAL_NCBI_CLEANASN
#include <accpubseq.h>
#endif
#define NLM_GENERATED_CODE_PROTO
#include <objmacro.h>
#include <macroapi.h>

#define CLEANASN_APP_VER "5.3"

CharPtr CLEANASN_APPLICATION = CLEANASN_APP_VER;

typedef struct sums {
  Int4          nucs;
  Int4          prts;
  Int4          recs;
} SumData, PNTR SumDataPtr;

typedef struct dbsums {
  SumData      genbank;
  SumData      embl;
  SumData      ddbj;
  SumData      refseq;
  SumData      other;
} DbSumData, PNTR DbSumPtr;

typedef struct counts {
  Int4          auth;
  Int4          bsec;
  Int4          clnr;
  Int4          gbbk;
  Int4          modr;
  Int4          move;
  Int4          norm;
  Int4          nucs;
  Int4          okay;
  Int4          othr;
  Int4          pack;
  Int4          prts;
  Int4          publ;
  Int4          recs;
  Int4          sloc;
  Int4          sort;
  Int4          ssec;
  Int4          titl;
} CountData, PNTR CountDataPtr;

typedef struct cleanflags {
  Char          buf [64];
  Int4          gi;
  Int2          year;
  Boolean       stripSerial;
  Boolean       isRefSeq;
  Boolean       isEmblDdbj;
  Boolean       batch;
  Boolean       binary;
  Boolean       compressed;
  Int2          type;
  CharPtr       results;
  CharPtr       outfile;
  CharPtr       firstfile;
  CharPtr       lastfile;
  Boolean       foundfirst;
  Boolean       foundlast;
  CharPtr       sourcedb;
  CharPtr       report;
  CharPtr       selective;
  ModType       ffmode;
  CharPtr       ffdiff;
  CharPtr       asn2flat;
  CharPtr       asnval;
  CharPtr       clean;
  CharPtr       modernize;
  CharPtr       link;
  CharPtr       feat;
  CharPtr       desc;
  CharPtr       mods;
  ValNodePtr    action_list;
  Boolean       taxon;
  CharPtr       pub;
  Int4          unpubcount;
  CountData     rawcounts;
  CountData     cumcounts;
  DbSumData     dbsums;
  AsnModulePtr  amp;
  AsnTypePtr    atp_bss;
  AsnTypePtr    atp_bsss;
  AsnTypePtr    atp_se;
  AsnTypePtr    atp_bsc;
  AsnTypePtr    bssp_atp;
  BioseqSet     bss;
  FILE          *logfp;
} CleanFlagData, PNTR CleanFlagPtr;

#ifdef INTERNAL_NCBI_CLEANASN
static CharPtr smartfetchproc = "SmartBioseqFetch";

static CharPtr smartfetchcmd = NULL;

extern Pointer ReadFromSmart (CharPtr accn, Uint2Ptr datatype, Uint2Ptr entityID);
extern Pointer ReadFromSmart (CharPtr accn, Uint2Ptr datatype, Uint2Ptr entityID)

{
  Char     cmmd [256];
  Pointer  dataptr;
  FILE*    fp;
  Char     path [PATH_MAX];

  if (datatype != NULL) {
    *datatype = 0;
  }
  if (entityID != NULL) {
    *entityID = 0;
  }
  if (StringHasNoText (accn)) return NULL;

  if (smartfetchcmd == NULL) {
    if (GetAppParam ("SEQUIN", "SMART", "FETCHSCRIPT", NULL, cmmd, sizeof (cmmd))) {
      smartfetchcmd = StringSaveNoNull (cmmd);
    }
  }
  if (smartfetchcmd == NULL) return NULL;

  TmpNam (path);

#ifdef OS_UNIX
  sprintf (cmmd, "csh %s %s > %s 2> /dev/null", smartfetchcmd, accn, path);
  system (cmmd);
#endif
#ifdef OS_MSWIN
  sprintf (cmmd, "%s %s -o %s", smartfetchcmd, accn, path);
  system (cmmd);
#endif

  fp = FileOpen (path, "r");
  if (fp == NULL) {
    FileRemove (path);
    return NULL;
  }
  dataptr = ReadAsnFastaOrFlatFile (fp, datatype, entityID, FALSE, FALSE, TRUE, FALSE);
  FileClose (fp);
  FileRemove (path);
  return dataptr;
}

static Int2 LIBCALLBACK SmartBioseqFetchFunc (Pointer data)

{
  BioseqPtr         bsp;
  Char              cmmd [256];
  Pointer           dataptr;
  Uint2             datatype;
  Uint2             entityID;
  FILE*             fp;
  OMProcControlPtr  ompcp;
  ObjMgrProcPtr     ompp;
  Char              path [PATH_MAX];
  SeqEntryPtr       sep = NULL;
  SeqIdPtr          sip;
  TextSeqIdPtr      tsip;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL) return OM_MSG_RET_ERROR;
  ompp = ompcp->proc;
  if (ompp == NULL) return OM_MSG_RET_ERROR;
  sip = (SeqIdPtr) ompcp->input_data;
  if (sip == NULL) return OM_MSG_RET_ERROR;

  if (sip->choice != SEQID_GENBANK) return OM_MSG_RET_ERROR;
  tsip = (TextSeqIdPtr) sip->data.ptrvalue;
  if (tsip == NULL || StringHasNoText (tsip->accession)) return OM_MSG_RET_ERROR;

  if (smartfetchcmd == NULL) {
    if (GetAppParam ("SEQUIN", "SMART", "FETCHSCRIPT", NULL, cmmd, sizeof (cmmd))) {
      smartfetchcmd = StringSaveNoNull (cmmd);
    }
  }
  if (smartfetchcmd == NULL) return OM_MSG_RET_ERROR;

  TmpNam (path);

#ifdef OS_UNIX
  sprintf (cmmd, "csh %s %s > %s 2> /dev/null", smartfetchcmd, tsip->accession, path);
  system (cmmd);
#endif
#ifdef OS_MSWIN
  sprintf (cmmd, "%s %s -o %s", smartfetchcmd, tsip->accession, path);
  system (cmmd);
#endif

  fp = FileOpen (path, "r");
  if (fp == NULL) {
    FileRemove (path);
    return OM_MSG_RET_ERROR;
  }
  dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, &entityID, FALSE, FALSE, TRUE, FALSE);
  FileClose (fp);
  FileRemove (path);

  if (dataptr == NULL) return OM_MSG_RET_OK;

  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep == NULL) return OM_MSG_RET_ERROR;
  bsp = BioseqFindInSeqEntry (sip, sep);
  ompcp->output_data = (Pointer) bsp;
  ompcp->output_entityID = ObjMgrGetEntityIDForChoice (sep);
  return OM_MSG_RET_DONE;
}

static Boolean SmartFetchEnable (void)

{
  ObjMgrProcLoad (OMPROC_FETCH, smartfetchproc, smartfetchproc,
                  OBJ_SEQID, 0, OBJ_BIOSEQ, 0, NULL,
                  SmartBioseqFetchFunc, PROC_PRIORITY_DEFAULT);
  return TRUE;
}

static CharPtr tpasmartfetchproc = "TPASmartBioseqFetch";

static CharPtr tpasmartfetchcmd = NULL;

extern Pointer ReadFromTPASmart (CharPtr accn, Uint2Ptr datatype, Uint2Ptr entityID);
extern Pointer ReadFromTPASmart (CharPtr accn, Uint2Ptr datatype, Uint2Ptr entityID)

{
  Char     cmmd [256];
  Pointer  dataptr;
  FILE*    fp;
  Char     path [PATH_MAX];

  if (datatype != NULL) {
    *datatype = 0;
  }
  if (entityID != NULL) {
    *entityID = 0;
  }
  if (StringHasNoText (accn)) return NULL;

  if (tpasmartfetchcmd == NULL) {
    if (GetAppParam ("SEQUIN", "TPASMART", "FETCHSCRIPT", NULL, cmmd, sizeof (cmmd))) {
      tpasmartfetchcmd = StringSaveNoNull (cmmd);
    }
  }
  if (tpasmartfetchcmd == NULL) return NULL;

  TmpNam (path);

#ifdef OS_UNIX
  sprintf (cmmd, "csh %s %s > %s 2> /dev/null", tpasmartfetchcmd, accn, path);
  system (cmmd);
#endif
#ifdef OS_MSWIN
  sprintf (cmmd, "%s %s -o %s", tpasmartfetchcmd, accn, path);
  system (cmmd);
#endif

  fp = FileOpen (path, "r");
  if (fp == NULL) {
    FileRemove (path);
    return NULL;
  }
  dataptr = ReadAsnFastaOrFlatFile (fp, datatype, entityID, FALSE, FALSE, TRUE, FALSE);
  FileClose (fp);
  FileRemove (path);
  return dataptr;
}

static Int2 LIBCALLBACK TPASmartBioseqFetchFunc (Pointer data)

{
  BioseqPtr         bsp;
  Char              cmmd [256];
  Pointer           dataptr;
  Uint2             datatype;
  Uint2             entityID;
  FILE*             fp;
  OMProcControlPtr  ompcp;
  ObjMgrProcPtr     ompp;
  Char              path [PATH_MAX];
  SeqEntryPtr       sep = NULL;
  SeqIdPtr          sip;
  TextSeqIdPtr      tsip;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL) return OM_MSG_RET_ERROR;
  ompp = ompcp->proc;
  if (ompp == NULL) return OM_MSG_RET_ERROR;
  sip = (SeqIdPtr) ompcp->input_data;
  if (sip == NULL) return OM_MSG_RET_ERROR;

  if (sip->choice != SEQID_TPG) return OM_MSG_RET_ERROR;
  tsip = (TextSeqIdPtr) sip->data.ptrvalue;
  if (tsip == NULL || StringHasNoText (tsip->accession)) return OM_MSG_RET_ERROR;

  if (tpasmartfetchcmd == NULL) {
    if (GetAppParam ("SEQUIN", "TPASMART", "FETCHSCRIPT", NULL, cmmd, sizeof (cmmd))) {
      tpasmartfetchcmd = StringSaveNoNull (cmmd);
    }
  }
  if (tpasmartfetchcmd == NULL) return OM_MSG_RET_ERROR;

  TmpNam (path);

#ifdef OS_UNIX
  sprintf (cmmd, "csh %s %s > %s 2> /dev/null", tpasmartfetchcmd, tsip->accession, path);
  system (cmmd);
#endif
#ifdef OS_MSWIN
  sprintf (cmmd, "%s %s -o %s", tpasmartfetchcmd, tsip->accession, path);
  system (cmmd);
#endif

  fp = FileOpen (path, "r");
  if (fp == NULL) {
    FileRemove (path);
    return OM_MSG_RET_ERROR;
  }
  dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, &entityID, FALSE, FALSE, TRUE, FALSE);
  FileClose (fp);
  FileRemove (path);

  if (dataptr == NULL) return OM_MSG_RET_OK;

  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep == NULL) return OM_MSG_RET_ERROR;
  bsp = BioseqFindInSeqEntry (sip, sep);
  ompcp->output_data = (Pointer) bsp;
  ompcp->output_entityID = ObjMgrGetEntityIDForChoice (sep);
  return OM_MSG_RET_DONE;
}

static Boolean TPASmartFetchEnable (void)

{
  ObjMgrProcLoad (OMPROC_FETCH, tpasmartfetchproc, tpasmartfetchproc,
                  OBJ_SEQID, 0, OBJ_BIOSEQ, 0, NULL,
                  TPASmartBioseqFetchFunc, PROC_PRIORITY_DEFAULT);
  return TRUE;
}

static CharPtr hupfetchproc = "HUPBioseqFetch";

static CharPtr hupfetchcmd = NULL;

extern Pointer ReadFromHUP (CharPtr accn, Uint2Ptr datatype, Uint2Ptr entityID);
extern Pointer ReadFromHUP (CharPtr accn, Uint2Ptr datatype, Uint2Ptr entityID)

{
  Char     cmmd [256];
  Pointer  dataptr;
  FILE*    fp;
  Char     path [PATH_MAX];

  if (datatype != NULL) {
    *datatype = 0;
  }
  if (entityID != NULL) {
    *entityID = 0;
  }
  if (StringHasNoText (accn)) return NULL;

  if (hupfetchcmd == NULL) {
    if (GetAppParam ("SEQUIN", "HUP", "FETCHSCRIPT", NULL, cmmd, sizeof (cmmd))) {
      hupfetchcmd = StringSaveNoNull (cmmd);
    }
  }
  if (hupfetchcmd == NULL) return NULL;

  TmpNam (path);

#ifdef OS_UNIX
  sprintf (cmmd, "csh %s %s > %s 2> /dev/null", hupfetchcmd, accn, path);
  system (cmmd);
#endif
#ifdef OS_MSWIN
  sprintf (cmmd, "%s %s -o %s", hupfetchcmd, accn, path);
  system (cmmd);
#endif

  fp = FileOpen (path, "r");
  if (fp == NULL) {
    FileRemove (path);
    return NULL;
  }
  dataptr = ReadAsnFastaOrFlatFile (fp, datatype, entityID, FALSE, FALSE, TRUE, FALSE);
  FileClose (fp);
  FileRemove (path);
  return dataptr;
}

static Int2 LIBCALLBACK HUPBioseqFetchFunc (Pointer data)

{
  BioseqPtr         bsp;
  Char              cmmd [256];
  Pointer           dataptr;
  Uint2             datatype;
  Uint2             entityID;
  FILE*             fp;
  OMProcControlPtr  ompcp;
  ObjMgrProcPtr     ompp;
  Char              path [PATH_MAX];
  SeqEntryPtr       sep = NULL;
  SeqIdPtr          sip;
  TextSeqIdPtr      tsip;

  ompcp = (OMProcControlPtr) data;
  if (ompcp == NULL) return OM_MSG_RET_ERROR;
  ompp = ompcp->proc;
  if (ompp == NULL) return OM_MSG_RET_ERROR;
  sip = (SeqIdPtr) ompcp->input_data;
  if (sip == NULL) return OM_MSG_RET_ERROR;

  if (sip->choice != SEQID_GENBANK) return OM_MSG_RET_ERROR;
  tsip = (TextSeqIdPtr) sip->data.ptrvalue;
  if (tsip == NULL || StringHasNoText (tsip->accession)) return OM_MSG_RET_ERROR;

  if (hupfetchcmd == NULL) {
    if (GetAppParam ("SEQUIN", "HUP", "FETCHSCRIPT", NULL, cmmd, sizeof (cmmd))) {
      hupfetchcmd = StringSaveNoNull (cmmd);
    }
  }
  if (hupfetchcmd == NULL) return OM_MSG_RET_ERROR;

  TmpNam (path);

#ifdef OS_UNIX
  sprintf (cmmd, "csh %s %s > %s 2> /dev/null", hupfetchcmd, tsip->accession, path);
  system (cmmd);
#endif
#ifdef OS_MSWIN
  sprintf (cmmd, "%s %s -o %s", hupfetchcmd, tsip->accession, path);
  system (cmmd);
#endif

  fp = FileOpen (path, "r");
  if (fp == NULL) {
    FileRemove (path);
    return OM_MSG_RET_ERROR;
  }
  dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, &entityID, FALSE, FALSE, TRUE, FALSE);
  FileClose (fp);
  FileRemove (path);

  if (dataptr == NULL) return OM_MSG_RET_OK;

  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep == NULL) return OM_MSG_RET_ERROR;
  bsp = BioseqFindInSeqEntry (sip, sep);
  ompcp->output_data = (Pointer) bsp;
  ompcp->output_entityID = ObjMgrGetEntityIDForChoice (sep);
  return OM_MSG_RET_DONE;
}

static Boolean HUPFetchEnable (void)

{
  ObjMgrProcLoad (OMPROC_FETCH, hupfetchproc, hupfetchproc,
                  OBJ_SEQID, 0, OBJ_BIOSEQ, 0, NULL,
                  HUPBioseqFetchFunc, PROC_PRIORITY_DEFAULT);
  return TRUE;
}
#endif

static void RemoveFeatUser (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  if (sfp == NULL) return;
  if (sfp->ext != NULL) {
    sfp->ext = UserObjectFree (sfp->ext);
  }
}

static void RemoveFeatDbxref (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  DbtagPtr    dbt;
  ValNodePtr  next, vnp;

  if (sfp == NULL) return;
  for (vnp = sfp->dbxref; vnp != NULL; vnp = next) {
    next = vnp->next;
    dbt = (DbtagPtr) vnp->data.ptrvalue;
    DbtagFree (dbt);
    MemFree (vnp);
  }
  sfp->dbxref = NULL;
}

typedef struct dummysmfedata {
  Int4  max;
  Int4  num_at_max;
} DummySmfeData, PNTR DummySmfePtr;

static Boolean LIBCALLBACK CADummySMFEProc (
  SeqFeatPtr sfp,
  SeqMgrFeatContextPtr context
)

{
  DummySmfePtr  dsp;
  Int4          len;

  if (sfp == NULL || context == NULL) return TRUE;
  dsp = context->userdata;
  if (dsp == NULL) return TRUE;

  len = SeqLocLen (sfp->location);
  if (len < dsp->max) {
    dsp->max = len;
    dsp->num_at_max = 1;
  } else if (len == dsp->max) {
    (dsp->num_at_max)++;
  }

  return TRUE;
}

static void RemoveUnnecGeneXref (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  Int2                 count;
  SeqFeatXrefPtr       curr, next;
  DummySmfeData        dsd;
  SeqMgrFeatContext    fcontext;
  SeqFeatXrefPtr PNTR  last;
  GeneRefPtr           grp, grpx;
  SeqFeatPtr           sfpx;
  CharPtr              syn1, syn2;

  if (sfp == NULL || sfp->data.choice == SEQFEAT_GENE) return;
  grp = SeqMgrGetGeneXref (sfp);
  if (grp == NULL || SeqMgrGeneIsSuppressed (grp)) return;
  sfpx = SeqMgrGetOverlappingGene (sfp->location, &fcontext);
  if (sfpx == NULL || sfpx->data.choice != SEQFEAT_GENE) return;
  grpx = (GeneRefPtr) sfpx->data.value.ptrvalue;
  if (grpx == NULL) return;

  if (StringDoesHaveText (grp->locus_tag) && StringDoesHaveText (grpx->locus_tag)) {
    if (StringICmp (grp->locus_tag, grpx->locus_tag) != 0) return;
  } else if (StringDoesHaveText (grp->locus) && StringDoesHaveText (grpx->locus)) {
    if (StringICmp (grp->locus, grpx->locus) != 0) return;
  } else if (grp->syn != NULL && grpx->syn != NULL) {
    syn1 = (CharPtr) grp->syn->data.ptrvalue;
    syn2 = (CharPtr) grpx->syn->data.ptrvalue;
    if (StringDoesHaveText (syn1) && StringDoesHaveText (syn2)) {
      if (StringICmp (syn1, syn2) != 0) return;
    }
  }

  MemSet ((Pointer) &dsd, 0, sizeof (DummySmfeData));
  dsd.max = INT4_MAX;
  dsd.num_at_max = 0;
  count = SeqMgrGetAllOverlappingFeatures (sfp->location, FEATDEF_GENE,
                                           NULL, 0, LOCATION_SUBSET,
                                           (Pointer) &dsd, CADummySMFEProc);

  if (dsd.num_at_max < 2) {
    last = (SeqFeatXrefPtr PNTR) &(sfp->xref);
    curr = sfp->xref;
    while (curr != NULL) {
      next = curr->next;
      if (curr->data.choice == SEQFEAT_GENE) {
        *last = next;
        curr->next = NULL;
        SeqFeatXrefFree (curr);
      } else {
        last = &(curr->next);
      }
      curr = next;
    }
  }
}

static void MarkTitles (
  SeqDescrPtr sdp,
  Pointer userdata
)

{
  ObjValNodePtr  ovn;

  if (sdp == NULL || sdp->choice != Seq_descr_title) return;
  if (sdp->extended == 0) return;
  ovn = (ObjValNodePtr) sdp;
  ovn->idx.deleteme = TRUE;
}

static void DoAutoDef (
  SeqEntryPtr sep,
  Uint2 entityID
)

{
  ValNodePtr                    defline_clauses = NULL;
  DeflineFeatureRequestList     feature_requests;
  Int4                          index;
  ValNodePtr                    modifier_indices = NULL;
  ModifierItemLocalPtr          modList;
  OrganismDescriptionModifiers  odmp;
  SeqEntryPtr                   oldscope;

  if (sep == NULL) return;
  if (entityID < 1) return;

  modList = MemNew (NumDefLineModifiers () * sizeof (ModifierItemLocalData));
  if (modList == NULL) return;

  InitFeatureRequests (&feature_requests);

  SetRequiredModifiers (modList);
  CountModifiers (modList, sep);

  InitOrganismDescriptionModifiers (&odmp, sep);

  RemoveNucProtSetTitles (sep);
  oldscope = SeqEntrySetScope (sep);

  BuildDefLineFeatClauseList (sep, entityID, &feature_requests,
                              DEFAULT_ORGANELLE_CLAUSE, FALSE, FALSE,
                              &defline_clauses);
  if (AreFeatureClausesUnique (defline_clauses)) {
    modifier_indices = GetModifierIndicesFromModList (modList);
  } else {
    modifier_indices = FindBestModifiers (sep, modList);
  }

  BuildDefinitionLinesFromFeatureClauseLists (defline_clauses, modList,
                                              modifier_indices, &odmp);
  DefLineFeatClauseListFree (defline_clauses);
  if (modList != NULL) {
    for (index = 0; index < NumDefLineModifiers (); index++) {
      ValNodeFree (modList [index].values_seen);
    }
    MemFree (modList);
  }
  modifier_indices = ValNodeFree (modifier_indices);

  ClearProteinTitlesInNucProts (entityID, NULL);
  InstantiateProteinTitles (entityID, NULL);
  /*
  RemovePopsetTitles (sep);
  */
  AddPopsetTitles (sep, &feature_requests, DEFAULT_ORGANELLE_CLAUSE, FALSE, FALSE);

  SeqEntrySetScope (oldscope);
}

static void LookupPubdesc (
  PubdescPtr pdp,
  Pointer userdata
)

{
  CitArtPtr        cap;
  MedlineEntryPtr  mep;
  PubmedEntryPtr   pep;
  Int4             pmid = 0;
  ValNodePtr       vnp;

  if (pdp == NULL) return;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    switch (vnp->choice) {
      case PUB_Muid :
        /* ignore obsolete muids */
        break;
      case PUB_PMid :
        pmid = vnp->data.intvalue;
        break;
      default :
        /* return on real pub */
        return;
        break;
    }
  }

  if (pmid == 0) return;

  pep = GetPubMedForUid (pmid);
  if (pep == NULL) return;
  mep = (MedlineEntryPtr) pep->medent;
  if (mep != NULL && mep->cit != NULL) {
    cap = AsnIoMemCopy ((Pointer) mep->cit,
                        (AsnReadFunc) CitArtAsnRead,
                        (AsnWriteFunc) CitArtAsnWrite);
    ValNodeAddPointer (&(pdp->pub), PUB_Article, (Pointer) cap);
  }

  PubmedEntryFree (pep);
}

static AuthListPtr AppendConsortiumToAuthList (
  AuthListPtr alp,
  CharPtr consortium
)

{
  AuthorPtr    ap;
  ValNodePtr   names;
  PersonIdPtr  pid;

  if (StringHasNoText (consortium)) return alp;
  if (alp == NULL) {
    alp = AuthListNew ();
    alp->choice = 1;
  }
  pid = PersonIdNew ();
  if (pid == NULL) return NULL;
  pid->choice = 5;
  pid->data = StringSave (consortium);
  ap = AuthorNew ();
  if (ap == NULL) return NULL;
  ap->name = pid;
  names = ValNodeAdd (&(alp->names));
  names->choice = 1;
  names->data.ptrvalue = ap;
  return alp;
}

static void ReplaceUnpub (
  PubdescPtr pdp,
  Pointer userdata
)

{
  CharPtr          authors = NULL, consortium = NULL;
  CitArtPtr        cap;
  CitGenPtr        cgp = NULL;
  Int2             count = 0;
  Boolean          hasUnpublished = FALSE;
  MedlineEntryPtr  mep;
  PubmedEntryPtr   pep;
  Int4             pmid;
  Int4Ptr          pmP;
  ValNodePtr       vnp = NULL;

  if (pdp == NULL || userdata == NULL) return;
  pmP = (Int4Ptr) userdata;
  pmid = *pmP;
  if (pmid == 0) return;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    count++;
    if (vnp->choice == PUB_Gen) {
      cgp = (CitGenPtr) vnp->data.ptrvalue;
      if (cgp != NULL) {
        if (StringICmp (cgp->cit, "Unpublished") == 0) {
          if (StringICmp (cgp->title, "Direct Submission") != 0) {
            hasUnpublished = TRUE;
          }
        }
      }
    } else if (vnp->choice == PUB_Muid || vnp->choice == PUB_PMid) {
      return;
    } else if (vnp->choice == PUB_Article || vnp->choice == PUB_Book || vnp->choice == PUB_Man) {
      return;
    }
  }

  if (! hasUnpublished) return;

  /* assume only a cit-gen for now */
  if (count != 1 || cgp == NULL) return;

  authors = GetAuthorsString (GENBANK_FMT, cgp->authors, &consortium, NULL, NULL);
  MemFree (authors);

  pep = GetPubMedForUid (pmid);
  if (pep == NULL) return;
  mep = (MedlineEntryPtr) pep->medent;
  if (mep != NULL && mep->cit != NULL) {
    /* remove cit-gen */
    CitGenFree (cgp);
    ValNodeFree (pdp->pub);
    pdp->pub = NULL;

    ValNodeAddInt (&(pdp->pub), PUB_PMid, pmid);
    cap = AsnIoMemCopy ((Pointer) mep->cit,
                        (AsnReadFunc) CitArtAsnRead,
                        (AsnWriteFunc) CitArtAsnWrite);
    ValNodeAddPointer (&(pdp->pub), PUB_Article, (Pointer) cap);
    if (consortium != NULL) {
      AppendConsortiumToAuthList (cap->authors, consortium);
    }
  }

  MemFree (consortium);

  PubmedEntryFree (pep);
}

static void CleanupLocation (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  BioseqPtr  bsp;
  SeqIntPtr  sintp;
  SeqLocPtr  slp;

  if (sfp == NULL || sfp->location == NULL) return;

  CleanUpSeqLoc (sfp->location);

  if (sfp->data.choice == SEQFEAT_REGION ||
      sfp->data.choice == SEQFEAT_SITE ||
      sfp->data.choice == SEQFEAT_BOND ||
      sfp->data.choice == SEQFEAT_PROT) {
    bsp = BioseqFind (SeqLocId (sfp->location));
    if (bsp != NULL && ISA_aa (bsp->mol)) {
      slp = SeqLocFindNext (sfp->location, NULL);
      while (slp != NULL) {
        if (slp->choice == SEQLOC_INT) {
          sintp = (SeqIntPtr) slp->data.ptrvalue;
          if (sintp != NULL) {
            if (sintp->strand != Seq_strand_unknown) {
              sintp->strand = Seq_strand_unknown;
            }
          }
        }
        slp = SeqLocFindNext (sfp->location, slp);
      }
    }
  }
}

static void CleanupMostRNAs (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  RnaRefPtr  rrp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA) return;
  rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
  if (rrp == NULL || rrp->type == 255) return;

  CleanUpSeqFeat (sfp, FALSE, FALSE, TRUE, FALSE, NULL);
}

static void CleanupRemainingRNAs (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  if (sfp == NULL || sfp->data.choice != SEQFEAT_RNA) return;

  CleanUpSeqFeat (sfp, FALSE, FALSE, TRUE, FALSE, NULL);
}

static void CleanupPubAuthors (
  PubdescPtr pdp,
  Pointer userdata
)

{
  if (pdp == NULL) return;

  CleanUpPubdescAuthors (pdp);
}

static void CleanupPubBody (
  PubdescPtr pdp,
  Pointer userdata
)

{
  CleanFlagPtr  cfp;

  if (pdp == NULL) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp == NULL) return;

  CleanUpPubdescBody (pdp, cfp->stripSerial);
}

static void ModGenes (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  ModernizeGeneFields (sfp);
}

static void ModRNAs (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  ModernizeRNAFields (sfp);
}

static void ModPCRs (
  BioSourcePtr biop,
  Pointer userdata
)

{
  ModernizePCRPrimers (biop);
}

static void FixBspMol (
  BioseqPtr bsp,
  Pointer userdata
)

{
  SeqDescrPtr  desc;
  MolInfoPtr   mip;

  if (bsp == NULL) return;

  desc = GetNextDescriptorUnindexed (bsp, Seq_descr_molinfo, NULL);
  if (desc == NULL || desc->choice != Seq_descr_molinfo) return;
  mip = (MolInfoPtr) desc->data.ptrvalue;
  if (mip == NULL) return;
  if (bsp->mol == 0) {
    switch (mip->biomol) {
      case MOLECULE_TYPE_GENOMIC :
        bsp->mol = Seq_mol_na;
        break;
      case MOLECULE_TYPE_PRE_MRNA :
      case MOLECULE_TYPE_MRNA :
      case MOLECULE_TYPE_RRNA :
      case MOLECULE_TYPE_TRNA :
      case MOLECULE_TYPE_SNRNA :
      case MOLECULE_TYPE_SCRNA :
      case MOLECULE_TYPE_CRNA :
      case MOLECULE_TYPE_SNORNA :
      case MOLECULE_TYPE_TRANSCRIBED_RNA :
      case MOLECULE_TYPE_NCRNA :
      case MOLECULE_TYPE_TMRNA :
        bsp->mol = Seq_mol_rna;
        break;
      case MOLECULE_TYPE_PEPTIDE :
        bsp->mol = Seq_mol_aa;
        break;
      case MOLECULE_TYPE_OTHER_GENETIC_MATERIAL :
        bsp->mol = Seq_mol_other;
        break;
      case MOLECULE_TYPE_GENOMIC_MRNA_MIX :
        bsp->mol = Seq_mol_na;
        break;
      default :
        break;
    }
  } else if (bsp->mol != Seq_mol_rna
             && (mip->biomol == MOLECULE_TYPE_CRNA || mip->biomol == MOLECULE_TYPE_MRNA)) {
    bsp->mol = Seq_mol_rna;
  }
}

static ByteStorePtr Se2Bs (
  SeqEntryPtr sep
)

{
  AsnIoBSPtr    aibp;
  ByteStorePtr  bs;

  if (sep == NULL) return NULL;

  bs = BSNew (1000);
  if (bs == NULL) return NULL;
  aibp = AsnIoBSOpen ("wb", bs);
  if (aibp == NULL || aibp->aip == NULL) return NULL;

  SeqEntryAsnWrite (sep, aibp->aip, NULL);

  AsnIoFlush (aibp->aip);
  AsnIoBSClose (aibp);

  return bs;
}

static ByteStorePtr Se2BsX (
  SeqEntryPtr sep
)

{
  AsnIoBSPtr    aibp;
  ByteStorePtr  bs;

  if (sep == NULL) return NULL;

  bs = BSNew (1000);
  if (bs == NULL) return NULL;
  aibp = AsnIoBSOpen ("wb", bs);
  if (aibp == NULL || aibp->aip == NULL) return NULL;

  aibp->aip->asn_no_newline = TRUE;
  aibp->aip->asn_alt_struct = TRUE;

  SeqEntryAsnWrite (sep, aibp->aip, NULL);

  AsnIoFlush (aibp->aip);
  AsnIoBSClose (aibp);

  return bs;
}

/*
static CharPtr Se2Str (
  SeqEntryPtr sep
)

{
  AsnIoBSPtr    aibp;
  ByteStorePtr  bs;
  CharPtr       str;

  if (sep == NULL) return NULL;

  bs = BSNew (1000);
  if (bs == NULL) return NULL;
  aibp = AsnIoBSOpen ("w", bs);
  if (aibp == NULL) return NULL;

  SeqEntryAsnWrite (sep, aibp->aip, NULL);

  AsnIoFlush (aibp->aip);
  AsnIoBSClose (aibp);

  str = BSMerge (bs, NULL);
  BSFree (bs);

  return str;
}
*/

typedef struct chgdata {
  Boolean       isRefSeq;
  Boolean       sgml;
  Boolean       cdscodon;
  Boolean       rubisco;
  Boolean       rbc;
  Boolean       its;
  Boolean       rnaother;
  Boolean       trnanote;
  Boolean       oldbiomol;
  Boolean       oldgbqual;
  Boolean       badDbxref;
  Boolean       refDbxref;
  Boolean       srcDbxref;
  Boolean       capDbxref;
  Boolean       oldDbxref;
  Boolean       privDbxref;
  Boolean       multDbxref;
  Boolean       rareDbxref;
  Boolean       badOrg;
  Boolean       rpt_unit_seq;
  Boolean       hasUnpublished;
  Boolean       hasPublished;
  Int4          protdesc;
  Int4          sfpnote;
  Int4          gbsource;
  Int4          cdsconf;
} ChangeData, PNTR ChangeDataPtr;

static Boolean IsRubisco (
  CharPtr name
)

{
  return (StringICmp (name, "rubisco large subunit") == 0 ||
          StringICmp (name, "rubisco small subunit") == 0);
}

static Boolean IsRbc (
  CharPtr name
)

{
  return (StringICmp (name, "RbcL") == 0 ||
          StringICmp (name, "RbcS") == 0);
}

static Boolean IsITS (
  CharPtr name
)

{
  return (StringICmp (name, "its1") == 0 ||
          StringICmp (name, "its 1") == 0 ||
          StringICmp (name, "its2") == 0 ||
          StringICmp (name, "its 2") == 0 ||
          StringICmp (name, "its3") == 0 ||
          StringICmp (name, "its 3") == 0 ||
          StringICmp (name, "Ribosomal DNA internal transcribed spacer 1") == 0 ||
          StringICmp (name, "Ribosomal DNA internal transcribed spacer 2") == 0 ||
          StringICmp (name, "Ribosomal DNA internal transcribed spacer 3") == 0 ||
          StringICmp (name, "internal transcribed spacer 1 (ITS1)") == 0 ||
          StringICmp (name, "internal transcribed spacer 2 (ITS2)") == 0 ||
          StringICmp (name, "internal transcribed spacer 3 (ITS3)") == 0);
}

static Boolean HasSgml (
  CharPtr str
)

{
  Int2  ascii_len;
  Char  buf [1024];

  if (StringHasNoText (str)) return FALSE;

  ascii_len = Sgml2AsciiLen (str);
  if (ascii_len + 2 > sizeof (buf)) return FALSE;

  Sgml2Ascii (str, buf, ascii_len + 1);
  if (StringCmp (str, buf) != 0) {
    return TRUE;
  }

  return FALSE;
}

static void LookForBadDbxref (
  ValNodePtr list,
  ChangeDataPtr cdp,
  Boolean isSource
)

{
  Boolean      cap;
  DbtagPtr     dp;
  CharPtr      good;
  ObjectIdPtr  oip;
  Boolean      ref;
  Boolean      src;
  CharPtr      str;
  ValNodePtr   vnp;

  if (list == NULL || cdp == NULL) return;

  for (vnp = list; vnp != NULL; vnp = vnp->next) {
    dp = (DbtagPtr) vnp->data.ptrvalue;
    if (dp != NULL && StringDoesHaveText (dp->db)) {

      oip = dp->tag;
      if (oip != NULL && StringDoesHaveText (oip->str)) {
        if (StringChr (oip->str, ':') != NULL) {
          cdp->multDbxref = TRUE;
        }
      }

      str = dp->db;
      if (StringICmp (str, "PID") == 0 ||
          StringICmp (str, "PIDg") == 0 ||
          StringICmp (str, "PIDd") == 0 ||
          StringICmp (str, "PIDe") == 0 ||
          StringICmp (str, "NID") == 0 ||
          StringICmp (str, "GI") == 0) {
        cdp->privDbxref = TRUE;
        continue;
      }
      if (StringICmp (str, "SWISS-PROT") == 0 ||
          StringICmp (str, "SWISSPROT") == 0 ||
          StringICmp (str, "SPTREMBL") == 0 ||
          StringICmp (str, "SUBTILIS") == 0 ||
          StringICmp (str, "MGD") == 0 ||
          StringCmp (str, "cdd") == 0 ||
          StringICmp (str, "TrEMBL") == 0 ||
          StringICmp (str, "LocusID") == 0 ||
          StringICmp (str, "MaizeDB") == 0 ||
          StringICmp (str, "UniProt/Swiss-Prot") == 0 ||
          StringICmp (str, "UniProt/TrEMBL") == 0 ||
          StringICmp (str, "Genew") == 0 ||
          StringICmp (str, "GENEDB") == 0 ||
          StringICmp (str, "GreengenesID") == 0 ||
          StringICmp (str, "HMPID") == 0 ||
          StringICmp (str, "IFO") == 0 ||
          StringICmp (str, "BHB") == 0 ||
          StringICmp (str, "BioHealthBase") == 0) {
        cdp->oldDbxref = TRUE;
        continue;
      }
      if (StringICmp (str, "ATCC(dna)") == 0 ||
          StringICmp (str, "ATCC(in host)") == 0 ||
          StringICmp (str, "BDGP_EST") == 0 ||
          StringICmp (str, "BDGP_INS") == 0 ||
          StringICmp (str, "CGNC") == 0 ||
          StringICmp (str, "CloneID") == 0 ||
          StringICmp (str, "ENSEMBL") == 0 ||
          StringICmp (str, "ESTLIB") == 0 ||
          StringICmp (str, "GDB") == 0 ||
          /*
          StringICmp (str, "GOA") == 0 ||
          */
          StringICmp (str, "IMGT/HLA") == 0 ||
          StringICmp (str, "PIR") == 0 ||
          StringICmp (str, "PSEUDO") == 0 ||
          StringICmp (str, "RZPD") == 0 ||
          StringICmp (str, "SoyBase") == 0 ||
          StringICmp (str, "UNILIB") == 0) {
        cdp->rareDbxref = TRUE;
        continue;
      }
      if (StringICmp (str, "MGD") == 0 || StringICmp (str, "MGI") == 0) {
        oip = dp->tag;
        if (oip != NULL && StringDoesHaveText (oip->str)) {
          str = oip->str;
          if (StringNICmp (str, "MGI:", 4) == 0 || StringNICmp (str, "MGD:", 4) == 0) {
            cdp->oldDbxref = TRUE;
            continue;
          }
        }
      } else if (StringICmp (str, "HPRD") == 0) {
        oip = dp->tag;
        if (oip != NULL && StringDoesHaveText (oip->str)) {
          str = oip->str;
          if (StringNICmp (str, "HPRD_", 5) == 0) {
            cdp->oldDbxref = TRUE;
            continue;
          }
        }
      }

      if (isSource && StringCmp (str, "taxon") == 0) continue;

      if (DbxrefIsValid (str, &ref, &src, &cap, &good)) {
        if (ref && (! cdp->isRefSeq)) {
          cdp->refDbxref = TRUE;
        }
        if (isSource && (! src)) {
          cdp->srcDbxref = TRUE;
        }
        if (cap) {
          cdp->capDbxref = TRUE;
        }
      } else {
        cdp->badDbxref = TRUE;
      }
    }
  }
}

static void ScoreFeature (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  BioSourcePtr   biop;
  ChangeDataPtr  cdp;
  Char           ch;
  CharPtr        comment;
  CdRegionPtr    crp;
  CharPtr        desc;
  GBQualPtr      gbq;
  GeneRefPtr     grp;
  CharPtr        name;
  OrgRefPtr      orp;
  ProtRefPtr     prp;
  CharPtr        ptr;
  Uint1          residue;
  RnaRefPtr      rrp;
  CharPtr        str;
  ValNodePtr     vnp;

  if (sfp == NULL) return;
   cdp = (ChangeDataPtr) userdata;
   if (cdp == NULL) return;

  comment = sfp->comment;
  if (StringDoesHaveText (comment)) {
    (cdp->sfpnote)++;
  }

  for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
    if (StringCmp (gbq->qual, "partial") == 0 ||
        StringCmp (gbq->qual, "evidence") == 0 ||
        StringCmp (gbq->qual, "exception") == 0 ||
        StringCmp (gbq->qual, "note") == 0 ||
        StringCmp (gbq->qual, "notes") == 0 ||
        StringCmp (gbq->qual, "comment") == 0 ||
        StringCmp (gbq->qual, "db_xref") == 0 ||
        StringCmp (gbq->qual, "gdb_xref") == 0 ||
        StringCmp (gbq->qual, "rpt_unit") == 0 ||
        StringCmp (gbq->qual, "pseudo") == 0 ||
        StringCmp (gbq->qual, "gene") == 0 ||
        StringCmp (gbq->qual, "codon_start") == 0 ||
        StringCmp (gbq->qual, "transposon") == 0 ||
        StringCmp (gbq->qual, "insertion_seq") == 0) {
      cdp->oldgbqual = TRUE;
    } else if (StringICmp (gbq->qual, "rpt_unit_seq") == 0) {
      if (StringHasNoText (gbq->val)) continue;
      ptr = gbq->val;
      ch = *ptr;
      while (ch != '\0') {
        if (IS_UPPER (ch)) {
          cdp->rpt_unit_seq = TRUE;
        }
        ptr++;
        ch = *ptr;
      }
    }
  }

  LookForBadDbxref (sfp->dbxref, cdp, FALSE);

  /* skip feature types that do not use data.value.ptrvalue */
  switch (sfp->data.choice) {
    case SEQFEAT_COMMENT:
    case SEQFEAT_BOND:
    case SEQFEAT_SITE:
    case SEQFEAT_PSEC_STR:
      return;
    default:
      break;
  }

  if (sfp->data.value.ptrvalue == NULL) return;

  switch (sfp->data.choice) {
    case SEQFEAT_GENE:
      grp = (GeneRefPtr) sfp->data.value.ptrvalue;
      if (HasSgml (grp->locus)) {
        cdp->sgml = TRUE;
      }
      if (HasSgml (grp->desc)) {
        cdp->sgml = TRUE;
      }
      for (vnp = grp->syn; vnp != NULL; vnp = vnp->next) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (StringHasNoText (str)) continue;
        if (HasSgml (str)) {
          cdp->sgml = TRUE;
        }
      }
      for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
        if (StringCmp (gbq->qual, "map") == 0 ||
            StringCmp (gbq->qual, "allele") == 0 ||
            StringCmp (gbq->qual, "locus_tag") == 0 ||
            StringCmp (gbq->qual, "old_locus_tag") == 0) {
          cdp->oldgbqual = TRUE;
        }
      }
      LookForBadDbxref (grp->db, cdp, FALSE);
      break;
    case SEQFEAT_CDREGION:
      crp = (CdRegionPtr) sfp->data.value.ptrvalue;
      if (crp->conflict) {
        (cdp->cdsconf)++;
      }
      for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
        if (StringCmp (gbq->qual, "codon") != 0) continue;
        if (StringHasNoText (gbq->val)) continue;
        cdp->cdscodon = TRUE;
      }
      for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
        if (StringCmp (gbq->qual, "product") == 0 ||
            StringCmp (gbq->qual, "function") == 0 ||
            StringCmp (gbq->qual, "EC_number") == 0 ||
            StringCmp (gbq->qual, "prot_note") == 0) {
          cdp->oldgbqual = TRUE;
        }
      }
      break;
    case SEQFEAT_PROT:
      prp = (ProtRefPtr) sfp->data.value.ptrvalue;
      desc = prp->desc;
      if (StringDoesHaveText (desc)) {
        (cdp->protdesc)++;
      }
      for (vnp = prp->name; vnp != NULL; vnp = vnp->next) {
        str = (CharPtr) vnp->data.ptrvalue;
        if (StringHasNoText (str)) continue;
        if (IsRubisco (str)) {
          cdp->rubisco = TRUE;
        }
        if (IsRbc (str)) {
          cdp->rbc = TRUE;
        }
      }
      for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
        if (StringCmp (gbq->qual, "product") == 0 ||
            StringCmp (gbq->qual, "function") == 0 ||
            StringCmp (gbq->qual, "EC_number") == 0 ||
            StringCmp (gbq->qual, "label") == 0 ||
            StringCmp (gbq->qual, "allele") == 0) {
          cdp->oldgbqual = TRUE;
        }
        if (StringCmp (gbq->qual, "standard_name") == 0 && prp->name == NULL) {
          cdp->oldgbqual = TRUE;
        }
      }
      LookForBadDbxref (prp->db, cdp, FALSE);
      break;
    case SEQFEAT_RNA :
      rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
      if (rrp->type == 255 && rrp->ext.choice == 1) {
        name = (CharPtr) rrp->ext.value.ptrvalue;
        if (StringCmp (name, "misc_RNA") == 0) {
          for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
            if (StringCmp (gbq->qual, "product") != 0) continue;
            name = gbq->val;
            if (StringHasNoText (name)) continue;
            if (IsITS (name)) {
              cdp->its = TRUE;
            }
          }
        } else if (StringCmp (name, "ncRNA") == 0 || StringCmp (name, "tmRNA") == 0) {
        } else {
          cdp->rnaother = TRUE;
          if (IsITS (name)) {
            cdp->its = TRUE;
          }
        }
      } else if (rrp->type == 3 && rrp->ext.choice == 2) {
        if (StringDoesHaveText (comment)) {
          if (StringNCmp (comment, "aa: ", 4) == 0) {
            comment += 4;
          }
          residue = FindTrnaAA3 (comment);
          if (residue > 0 && residue != 255) {
            cdp->trnanote = TRUE;
          }
          residue = FindTrnaAA (comment);
          if (residue > 0 && residue != 255) {
            cdp->trnanote = TRUE;
          }
        }
      }
      for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
        if (StringCmp (gbq->qual, "product") == 0 ||
            StringCmp (gbq->qual, "ncRNA_class") == 0 ||
            StringCmp (gbq->qual, "tag_peptide") == 0) {
          cdp->oldgbqual = TRUE;
        }
      }
      break;
    case SEQFEAT_ORG :
      orp = (OrgRefPtr) sfp->data.value.ptrvalue;
      LookForBadDbxref (orp->db, cdp, TRUE);
      cdp->badOrg = TRUE;
      break;
    case SEQFEAT_BIOSRC :
      biop = (BioSourcePtr) sfp->data.value.ptrvalue;
      orp = biop->org;
      if (orp != NULL) {
        LookForBadDbxref (orp->db, cdp, TRUE);
      }
    default:
      break;
  }
}

static void ScoreDescriptor (
  SeqDescrPtr sdp,
  Pointer userdata
)

{
  BioSourcePtr   biop;
  ChangeDataPtr  cdp;
  GBBlockPtr     gbp;
  MolInfoPtr     mip;
  OrgRefPtr      orp;

  if (sdp == NULL) return;
  cdp = (ChangeDataPtr) userdata;
  if (cdp == NULL) return;

  switch (sdp->choice) {
    case Seq_descr_genbank :
      gbp = (GBBlockPtr) sdp->data.ptrvalue;
      if (gbp != NULL) {
        if (StringDoesHaveText (gbp->source)) {
          (cdp->gbsource)++;
        }
      }
      break;
    case Seq_descr_molinfo :
      mip = (MolInfoPtr) sdp->data.ptrvalue;
      if (mip != NULL) {
        switch (mip->biomol) {
          case MOLECULE_TYPE_SNRNA:
          case MOLECULE_TYPE_SCRNA:
          case MOLECULE_TYPE_SNORNA:
            cdp->oldbiomol = TRUE;
            break;
          default :
            break;
        }
      }
      break;
    case Seq_descr_org :
      orp = (OrgRefPtr) sdp->data.ptrvalue;
      if (orp != NULL) {
        LookForBadDbxref (orp->db, cdp, TRUE);
      }
      cdp->badOrg = TRUE;
      break;
    case Seq_descr_source :
      biop = (BioSourcePtr) sdp->data.ptrvalue;
      if (biop != NULL) {
        orp = biop->org;
        if (orp != NULL) {
          LookForBadDbxref (orp->db, cdp, TRUE);
        }
      }
      break;
    default :
      break;
  }
}

static void CheckForUnpubPub (
  PubdescPtr pdp,
  Pointer userdata
)

{
  ChangeDataPtr  cdp;
  CitGenPtr      cgp;
  ValNodePtr     vnp;

  if (pdp == NULL) return;
  cdp = (ChangeDataPtr) userdata;
  if (cdp == NULL) return;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == PUB_Gen) {
      cgp = (CitGenPtr) vnp->data.ptrvalue;
      if (cgp != NULL) {
        if (StringICmp (cgp->cit, "Unpublished") == 0) {
          if (StringICmp (cgp->title, "Direct Submission") != 0) {
            cdp->hasUnpublished = TRUE;
          }
        }
      }
    } else if (vnp->choice == PUB_Muid || vnp->choice == PUB_PMid) {
      cdp->hasPublished = TRUE;
    } else if (vnp->choice == PUB_Article || vnp->choice == PUB_Book || vnp->choice == PUB_Man) {
      cdp->hasPublished = TRUE;
    }
  }
}

static void CheckForChanges (
  SeqEntryPtr sep,
  ChangeDataPtr cdp
)

{
  if (sep == NULL || cdp == NULL) return;

  VisitFeaturesInSep (sep, (Pointer) cdp, ScoreFeature);
  VisitDescriptorsInSep (sep, (Pointer) cdp, ScoreDescriptor);
  VisitPubdescsInSep (sep, (Pointer) cdp, CheckForUnpubPub);
}

static void StripBadProtTitles (
  BioseqPtr bsp,
  Pointer userdata
)

{
  CharPtr            buf;
  size_t             buflen = 1001;
  ObjValNodePtr      ovp;
  SeqIdPtr           sip;
  CharPtr            title;
  ValNodePtr         vnp;

  if (bsp == NULL) return;
  if (! ISA_aa (bsp->mol)) return;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_OTHER) return;
  }

  vnp = BioseqGetSeqDescr (bsp, Seq_descr_title, NULL);
  if (vnp == NULL) return;
  title = (CharPtr) vnp->data.ptrvalue;
  if (StringHasNoText (title)) return;

  buf = MemNew (sizeof (Char) * (buflen + 1));
  if (buf == NULL) return;

  if (NewCreateDefLineBuf (NULL, bsp, buf, buflen, TRUE, FALSE)) {
    if (StringICmp (buf, title) != 0) {
      if (vnp->extended != 0) {
        ovp = (ObjValNodePtr) vnp;
        ovp->idx.deleteme = TRUE;
      }
    }
  }

  MemFree (buf);
}

static void BadProtTitleProc (
  SeqEntryPtr sep,
  Pointer mydata,
  Int4 index,
  Int2 indent
)

{
  BioseqSetPtr  bssp;

  if (sep == NULL) return;
  if (! IS_Bioseq_set (sep)) return;
  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  if (bssp->_class != BioseqseqSet_class_nuc_prot) return;
  VisitBioseqsInSep (sep, NULL, StripBadProtTitles);
}

static void BSSaveToFile (
  ByteStorePtr bs,
  CharPtr path
)

{
  Byte  buf [256];
  Int4  count;
  FILE  *fp;

  if (bs == NULL || StringHasNoText (path)) return;

  fp = FileOpen (path, "w");
  if (fp != NULL) {
    Nlm_BSSeek (bs, 0, SEEK_SET);
    count = BSRead (bs, buf, sizeof (buf));
    while (count > 0) {
      FileWrite (buf, count, 1, fp);
      count = BSRead (bs, buf, sizeof (buf));
    }
    FileClose (fp);
  }
}

typedef struct diffblock {
  ValNodePtr  head;
  ValNodePtr  tail;
} DiffBlock, PNTR DiffBlockPtr;

static void RecordDiffBlock (
  DiffBlockPtr dbp,
  CharPtr str
)

{
  ValNodePtr  vnp;

  if (dbp == NULL || StringHasNoText (str)) return;

  vnp = ValNodeCopyStr (&(dbp->tail), 0, str);
  if (dbp->head == NULL) {
    dbp->head = vnp;
  }
  dbp->tail = vnp;
}

static void WriteDiffBlock (
  DiffBlockPtr dbp,
  FILE *fp
)

{
  Char        ch;
  Int2        idx;
  Int2        margin = INT2_MAX;
  CharPtr     ptr;
  Int2        spaces;
  CharPtr     str;
  ValNodePtr  vnp;

  if (dbp == NULL || dbp->head == NULL || fp == NULL) return;

  for (vnp = dbp->head; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) continue;
    ch = str [0];
    if (ch == '<' || ch == '>') {
      ptr = str + 1;
      ch = *ptr;
      spaces = 0;
      while (ch == ' ') {
        spaces++;
        ptr++;
        ch = *ptr;
      }
      if (spaces < margin) {
        margin = spaces;
      }
    }
  }

  if (margin > 80) {
    margin = 80;
  }

  for (vnp = dbp->head; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) continue;
    ch = str [0];
    if (ch == '<' || ch == '>') {
      ptr = str + 1;
      ch = *ptr;
      idx = 0;
      while (idx < margin && ch == ' ') {
        idx++;
        ptr++;
        ch = *ptr;
      }
      fprintf (fp, "%c %s\n", str [0], ptr);
    } else if (ch == '-') {
      fprintf (fp, "---\n");
    } else if (ch == '=') {
      fprintf (fp, "===\n");
    }
  }

  fprintf (fp, "\n");
  fflush (fp);
}

static void ResetDiffBlock (
  DiffBlockPtr dbp
)

{
  if (dbp == NULL) return;

  dbp->head = ValNodeFreeData (dbp->head);
  dbp->tail = NULL;
}

static void ReportAsnDiffs (
  FILE *logfp,
  CharPtr id,
  ByteStorePtr bs1,
  ByteStorePtr bs2)

{
#ifdef OS_UNIX
  Char       ch;
  Char       cmmd [512];
  DiffBlock  db;
  int        diff;
  FileCache  fc;
  FILE       *fp;
  Char       line [512];
  Char       path1 [PATH_MAX];
  Char       path2 [PATH_MAX];
  Char       path3 [PATH_MAX];
  CharPtr    str;

  if (logfp == NULL || StringHasNoText (id)) return;
  if (bs1 == NULL || bs2 == NULL) return;

  TmpNam (path1);
  TmpNam (path2);
  TmpNam (path3);

  BSSaveToFile (bs1, path1);
  BSSaveToFile (bs2, path2);

  db.head = NULL;
  db.tail = NULL;

  sprintf (cmmd, "diff -b -h %s %s > %s", path1, path2, path3);
  diff = system (cmmd);

  if (diff > 0) {
    fp = FileOpen (path3, "r");
    if (fp != NULL) {
      fprintf (logfp, "\n\n%s\n\n", id);
      fflush (logfp);
      if (FileCacheSetup (&fc, fp)) {
        str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
        while (str != NULL) {
          ch = line [0];
          if (ch == '<' || ch == '>') {
            RecordDiffBlock (&db, line);
          } else if (ch == '-') {
            RecordDiffBlock (&db, "---");
          } else if (IS_DIGIT (ch)) {
            WriteDiffBlock (&db, logfp);
            ResetDiffBlock (&db);
            RecordDiffBlock (&db, "===");
          } else if (StringHasNoText (str)) {
            WriteDiffBlock (&db, logfp);
            ResetDiffBlock (&db);
          }
          str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
        }
        WriteDiffBlock (&db, logfp);
        ResetDiffBlock (&db);
      }
      fprintf (logfp, "//\n\n");
      FileClose (fp);
    }
  }

  sprintf (cmmd, "rm %s; rm %s; rm %s", path1, path2, path3);
  system (cmmd);
#endif
}

static void LogPopSetTitleResults (FILE *fp, PopSetRetroStatPtr stat, SeqEntryPtr sep)
{
  ValNode vn;
  CharPtr txt;

  if (fp == NULL || stat == NULL || sep == NULL || !IS_Bioseq_set (sep)) {
    return;
  }
  /* if test wasn't run, return */
  if (stat->feature_clause == 0 && stat->common_title == 0 && stat->uncalculatable == 0) {
    return;
  }

  MemSet (&vn, 0, sizeof (ValNode));
  vn.choice = OBJ_BIOSEQSET;
  vn.data.ptrvalue = sep->data.ptrvalue;
  txt = GetDiscrepancyItemText (&vn);
  fprintf (fp, "AutoDefPopSetResults: %s", txt);
  txt = MemFree (txt);
  fprintf (fp, "\tFeature Clause: %d\n\tCommon Title: %d\n\tUncalculatable: %d\n",
           stat->feature_clause, stat->common_title, stat->uncalculatable);
}

static void DoAsnDiffReport (
  SeqEntryPtr sep,
  CleanFlagPtr cfp
)

{
  ByteStorePtr         bs = NULL, tmp = NULL;
  Uint2                entityID;
  Boolean              okay = FALSE;
  PopSetRetroStatData  stat;

  if (sep == NULL || cfp == NULL) return;

  MemSet ((Pointer) &stat, 0, sizeof (PopSetRetroStatData));

  RemoveAllNcbiCleanupUserObjects (sep);

  entityID = ObjMgrGetEntityIDForChoice (sep);

  /* Capital letters avoid unwanted diffs on issues already fixed in ID */

  if (StringChr (cfp->selective, 'S') != NULL) {
    SeriousSeqEntryCleanup (sep, NULL, NULL);
    RemoveAllNcbiCleanupUserObjects (sep);
    NormalizeDescriptorOrder (sep);
  }
  if (StringChr (cfp->selective, 'B') != NULL) {
    BasicSeqEntryCleanup (sep);
  }
  if (StringChr (cfp->selective, 'A') != NULL) {
    VisitPubdescsInSep (sep, NULL, CleanupPubAuthors);
  }
  if (StringChr (cfp->selective, 'P') != NULL) {
    VisitPubdescsInSep (sep, (Pointer) cfp, CleanupPubBody);
  }
  if (StringChr (cfp->selective, 'L') != NULL) {
    VisitFeaturesInSep (sep, NULL, CleanupLocation);
  }
  if (StringChr (cfp->selective, 'R') != NULL) {
    VisitFeaturesInSep (sep, NULL, CleanupMostRNAs);
    VisitFeaturesInSep (sep, NULL, CleanupRemainingRNAs);
    VisitFeaturesInSep (sep, NULL, ModRNAs);
  }
  if (StringChr (cfp->selective, 'Q') != NULL) {
    SortSeqEntryQualifiers (sep);
  }
  if (StringChr (cfp->selective, 'G') != NULL) {
    EntryChangeGBSource (sep);
    EntryCheckGBBlock (sep);
  }
  if (StringChr (cfp->selective, 'K') != NULL) {
    MoveFeatsFromPartsSet (sep);
    move_cds_ex (sep, TRUE);
  }
  if (StringChr (cfp->selective, 'M') != NULL) {
    SeqEntryPubsAsn4 (sep, cfp->isEmblDdbj);
  }
  if (StringChr (cfp->selective, 'O') != NULL) {
    SeqEntryPubsAsn4Ex (sep, cfp->isEmblDdbj, FALSE);
  }
  if (StringChr (cfp->selective, 'D') != NULL) {
    SeqMgrIndexFeatures (entityID, 0);
    DoAutoDef (sep, entityID);
  }
  if (StringChr (cfp->selective, 'E') != NULL) {
    SeqMgrIndexFeatures (entityID, 0);
    PopSetAutoDefRetro (sep, &stat);
    LogPopSetTitleResults (cfp->logfp, &stat, sep);
  }

  NormalizeDescriptorOrder (sep);

  /* Look for change in single issue */

  bs = Se2BsX (sep);
  if (StringHasNoText (cfp->selective)) {
    okay = TRUE;
  } else if (StringChr (cfp->selective, 's') != NULL) {
    SeriousSeqEntryCleanup (sep, NULL, NULL);
    RemoveAllNcbiCleanupUserObjects (sep);
    NormalizeDescriptorOrder (sep);
    tmp = Se2BsX (sep);
    if (! BSEqual (bs, tmp)) {
      okay = TRUE;
    }
    tmp = BSFree (tmp);
  } else if (StringChr (cfp->selective, 'b') != NULL) {
    BasicSeqEntryCleanup (sep);
    tmp = Se2BsX (sep);
    if (! BSEqual (bs, tmp)) {
      okay = TRUE;
    }
    tmp = BSFree (tmp);
  } else if (StringChr (cfp->selective, 'a') != NULL) {
    VisitPubdescsInSep (sep, NULL, CleanupPubAuthors);
    tmp = Se2BsX (sep);
    if (! BSEqual (bs, tmp)) {
      okay = TRUE;
    }
    tmp = BSFree (tmp);
  } else if (StringChr (cfp->selective, 'p') != NULL) {
    VisitPubdescsInSep (sep, (Pointer) cfp, CleanupPubBody);
    tmp = Se2BsX (sep);
    if (! BSEqual (bs, tmp)) {
      okay = TRUE;
    }
    tmp = BSFree (tmp);
  } else if (StringChr (cfp->selective, 'l') != NULL) {
    VisitFeaturesInSep (sep, NULL, CleanupLocation);
    tmp = Se2BsX (sep);
    if (! BSEqual (bs, tmp)) {
      okay = TRUE;
    }
    tmp = BSFree (tmp);
  } else if (StringChr (cfp->selective, 'r') != NULL) {
    VisitFeaturesInSep (sep, NULL, CleanupMostRNAs);
    VisitFeaturesInSep (sep, NULL, CleanupRemainingRNAs);
    VisitFeaturesInSep (sep, NULL, ModRNAs);
    tmp = Se2BsX (sep);
    if (! BSEqual (bs, tmp)) {
      okay = TRUE;
    }
    tmp = BSFree (tmp);
  } else if (StringChr (cfp->selective, 'q') != NULL) {
    SortSeqEntryQualifiers (sep);
    tmp = Se2BsX (sep);
    if (! BSEqual (bs, tmp)) {
      okay = TRUE;
    }
    tmp = BSFree (tmp);
  } else if (StringChr (cfp->selective, 'g') != NULL) {
    EntryChangeGBSource (sep);
    EntryCheckGBBlock (sep);
    tmp = Se2BsX (sep);
    if (! BSEqual (bs, tmp)) {
      okay = TRUE;
    }
    tmp = BSFree (tmp);
  } else if (StringChr (cfp->selective, 'k') != NULL) {
    MoveFeatsFromPartsSet (sep);
    move_cds_ex (sep, TRUE);
    tmp = Se2BsX (sep);
    if (! BSEqual (bs, tmp)) {
      okay = TRUE;
    }
    tmp = BSFree (tmp);
  } else if (StringChr (cfp->selective, 'm') != NULL) {
    SeqEntryPubsAsn4 (sep, cfp->isEmblDdbj);
    NormalizeDescriptorOrder (sep);
    tmp = Se2BsX (sep);
    if (! BSEqual (bs, tmp)) {
      okay = TRUE;
    }
    tmp = BSFree (tmp);
  } else if (StringChr (cfp->selective, 'o') != NULL) {
    SeqEntryPubsAsn4Ex (sep, cfp->isEmblDdbj, FALSE);
    NormalizeDescriptorOrder (sep);
    tmp = Se2BsX (sep);
    if (! BSEqual (bs, tmp)) {
      okay = TRUE;
    }
    tmp = BSFree (tmp);
  } else if (StringChr (cfp->selective, 'd') != NULL) {
    SeqMgrIndexFeatures (entityID, 0);
    DoAutoDef (sep, entityID);
    NormalizeDescriptorOrder (sep);
    tmp = Se2BsX (sep);
    if (! BSEqual (bs, tmp)) {
      okay = TRUE;
    }
    tmp = BSFree (tmp);
  } else if (StringChr (cfp->selective, 'e') != NULL) {
    SeqMgrIndexFeatures (entityID, 0);
    PopSetAutoDefRetro (sep, &stat);
    LogPopSetTitleResults (cfp->logfp, &stat, sep);
    NormalizeDescriptorOrder (sep);
    if (stat.title_added) {
      okay = TRUE;
    }
  } else {
    okay = TRUE;
  }

  /* Report incremental diff */

  if (okay) {
    SeriousSeqEntryCleanup (sep, NULL, NULL);
    RemoveAllNcbiCleanupUserObjects (sep);
    NormalizeDescriptorOrder (sep);
    tmp = Se2BsX (sep);
    if (! BSEqual (bs, tmp)) {
      if (cfp->logfp != NULL) {
        ReportAsnDiffs (cfp->logfp, cfp->buf, bs, tmp);
      }
    }
    tmp = BSFree (tmp);
  }

  BSFree (bs);
}

static void GetPopPhyMutEcoType (
  BioseqSetPtr bssp,
  Pointer userdata
)

{
  CharPtr PNTR  set_type;

  if (bssp == NULL || userdata == NULL) return;
  set_type = (CharPtr PNTR) userdata;

  if (bssp->_class >= BioseqseqSet_class_mut_set && bssp->_class <= BioseqseqSet_class_eco_set) {
    switch (bssp->_class) {
      case BioseqseqSet_class_mut_set :
        *set_type = "MUT";
        break;
      case BioseqseqSet_class_pop_set :
        *set_type = "POP";
        break;
      case BioseqseqSet_class_phy_set :
        *set_type = "PHY";
        break;
      case BioseqseqSet_class_eco_set :
        *set_type = "ECO";
        break;
      default :
        break;
    }
  }
}

/*
static void ChromosomeScanCallback (
  BioseqPtr bsp,
  Pointer userdata
)

{
  BioSourcePtr       biop = NULL;
  CleanFlagPtr       cfp;
  Boolean            complete = FALSE;
  SeqMgrDescContext  dcontext;
  CharPtr            genome = NULL;
  MolInfoPtr         mip = NULL;
  OrgNamePtr         onp;
  OrgRefPtr          orp;
  SeqDescrPtr        sdp;

  if (bsp == NULL || userdata == NULL) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp->logfp == NULL) return;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_source, &dcontext);
  if (sdp != NULL) {
    biop = (BioSourcePtr) sdp->data.ptrvalue;
  }
  if (biop == NULL) return;
  if (biop->genome == GENOME_chromosome) {
    genome = "CHROM";
  } else if (biop->genome == GENOME_genomic) {
    genome = "GENOM";
  } else if (biop->genome == GENOME_unknown) {
    genome = "UNKWN";
  } else {
    return;
  }
  orp = biop->org;
  if (orp == NULL) return;
  onp = orp->orgname;
  if (onp == NULL) return;
  if (StringNICmp (onp->lineage, "Eukaryota; ", 11) != 0) return;

  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
  if (sdp != NULL) {
    mip = (MolInfoPtr) sdp->data.ptrvalue;
  }
  if (mip != NULL && mip->completeness == 1) {
    complete = TRUE;
  }

  if (cfp->logfp != NULL) {
    if (complete) {
      fprintf (cfp->logfp, "COMP %s %s\n", genome, cfp->buf);
    } else {
      fprintf (cfp->logfp, "PART %s %s\n", genome, cfp->buf);
    }
    fflush (cfp->logfp);
  }
}
*/

/*
static void PopPhyMutEcoScanCallback (
  BioseqSetPtr bssp,
  Pointer userdata
)

{
  CleanFlagPtr  cfp;
  SeqDescrPtr   sdp;
  CharPtr       set_type = NULL;
  CharPtr       title = NULL;

  if (bssp == NULL || userdata == NULL) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp->logfp == NULL) return;

  if (bssp->_class >= BioseqseqSet_class_mut_set && bssp->_class <= BioseqseqSet_class_eco_set) {
    GetPopPhyMutEcoType (bssp, (Pointer) &set_type);
    switch (bssp->_class) {
      case BioseqseqSet_class_mut_set :
        set_type = "MUT";
        break;
      case BioseqseqSet_class_pop_set :
        set_type = "POP";
        break;
      case BioseqseqSet_class_phy_set :
        set_type = "PHY";
        break;
      case BioseqseqSet_class_eco_set :
        set_type = "ECO";
        break;
      default :
        break;
    }
    for (sdp = bssp->descr; sdp != NULL; sdp = sdp->next) {
      if (sdp->choice != Seq_descr_title) continue;
      title = (CharPtr) sdp->data.ptrvalue;
    }
    if (StringHasNoText (set_type)) {
      set_type = "UNK";
    }
    if (StringDoesHaveText (title)) {
      fprintf (cfp->logfp, "%s set TITLE %s\n", set_type, cfp->buf);
    } else {
      fprintf (cfp->logfp, "%s set ANNON %s\n", set_type, cfp->buf);
    }
  }
}
*/

/*
static void PseudoScanCallback (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  CleanFlagPtr  cfp;
  CharPtr       str;

  if (sfp == NULL || userdata == NULL) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp->logfp == NULL) return;

  if (! sfp->pseudo) return;

  str = FindKeyFromFeatDefType (sfp->idx.subtype, FALSE);
  if (StringHasNoText (str)) {
    str = "?";
  }

  fprintf (cfp->logfp, "%s Pseudo %s\n", cfp->buf, str);
}
*/

/*
static void PubScanCallback (
  PubdescPtr pdp,
  Pointer userdata
)

{
  CitArtPtr     cap;
  CleanFlagPtr  cfp;
  CitJourPtr    cjp;
  ImprintPtr    imp;
  ValNodePtr    vnp;

  if (pdp == NULL) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp == NULL || cfp->logfp == NULL) return;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == PUB_Article) {
      cap = (CitArtPtr) vnp->data.ptrvalue;
      if (cap != NULL && cap->from == 1) {
        cjp = (CitJourPtr) cap->fromptr;
        if (cjp != NULL) {
          imp = cjp->imp;
          if (imp != NULL) {
            if (imp->part_sup != NULL && StringHasNoText (imp->part_sup)) {
              fprintf (cfp->logfp, "EMPTY PART_SUP %s\n", cfp->buf);
            }
            if (StringDoesHaveText (imp->issue)) {
              if (StringNICmp (imp->issue, "PT ", 3) == 0) {
                fprintf (cfp->logfp, "IMP_PT %s\n", cfp->buf);
              }
              if (StringDoesHaveText (imp->part_sup)) {
                fprintf (cfp->logfp, "SUPN %s\n", cfp->buf);
              }
              if (StringDoesHaveText (imp->part_supi)) {
                fprintf (cfp->logfp, "SUPI %s\n", cfp->buf);
              }
            }
          }
        }
      }
    }
  }
}
*/

static void DoASNScan (
  SeqEntryPtr sep,
  CleanFlagPtr cfp
)

{
  /*
  ByteStorePtr  bs = NULL, tmp = NULL;

  if (sep == NULL || cfp == NULL || cfp->logfp == NULL) return;
  */

  /*
  VisitBioseqsInSep (sep, (Pointer) cfp, ChromosomeScanCallback);
  */

  /*
  VisitSetsInSep (sep, (Pointer) cfp, PopPhyMutEcoScanCallback);
  */

  /*
  SeriousSeqEntryCleanup (sep, NULL, NULL);
  RemoveAllNcbiCleanupUserObjects (sep);
  NormalizeDescriptorOrder (sep);

  bs = Se2Bs (sep);

  SeqEntryPubsAsn4Ex (sep, cfp->isEmblDdbj, TRUE);
  NormalizeDescriptorOrder (sep);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    fprintf (cfp->logfp, "%s\n", cfp->buf);
    fflush (cfp->logfp);
  }
  BSFree (bs);
  bs = tmp;

  BSFree (bs);
  */

  /*
  VisitFeaturesInSep (sep, (Pointer) cfp, PseudoScanCallback);
  */

  /*
  VisitPubdescsInSep (sep, (Pointer) cfp, PubScanCallback);
  */
}

static void DoASNReport (
  SeqEntryPtr sep,
  CleanFlagPtr cfp,
  Boolean dossec,
  Boolean quick,
  Boolean popphymutdef
)

{
  Boolean              auth = FALSE, bsec = FALSE, clnr = FALSE, gbbk = FALSE,
                       modr = FALSE, move = FALSE, norm = FALSE, othr = FALSE,
                       pack = FALSE, publ = FALSE, ssec = FALSE, sloc = FALSE,
                       sort = FALSE, titl = FALSE, popphymuttitle = FALSE,
                       chkseg = FALSE, chknps = FALSE, dnarna = FALSE,
                       ncbiusrobj = FALSE;
  ByteStorePtr         bs = NULL, tmp = NULL;
  ChangeData           cdbefore, cdafter;
  Int2                 chk_nuc_prot_val = 0;
  Uint2                entityID;
  CharPtr              set_type = NULL;
  PopSetRetroStatData  stat;

  if (sep == NULL || cfp == NULL) return;

  if (FindNcbiCleanupUserObject (sep) != NULL) {
    ncbiusrobj = TRUE;
  }

  MemSet ((Pointer) &stat, 0, sizeof (PopSetRetroStatData));

  RemoveAllNcbiCleanupUserObjects (sep);

  if (popphymutdef) {
    VisitSetsInSep (sep, (Pointer) &set_type, GetPopPhyMutEcoType);
    if (StringHasNoText (set_type)) return;

    entityID = ObjMgrGetEntityIDForChoice (sep);

    SeqMgrIndexFeatures (entityID, 0);
    PopSetAutoDefRetro (sep, &stat);
    LogPopSetTitleResults (cfp->logfp, &stat, sep);

    if (stat.title_added) {
      popphymuttitle = TRUE;
    }

    if (popphymuttitle) {
      if (cfp->logfp != NULL) {
        if (StringHasNoText (set_type)) {
          set_type = "UNK";
        }
        fprintf (cfp->logfp, "%s_SET %s\n", set_type, cfp->buf);
        fflush (cfp->logfp);
      }
    }

    return;
  }

  if (quick) {
    bs = Se2Bs (sep);

    NormalizeDescriptorOrder (sep);
    tmp = Se2Bs (sep);
    if (! BSEqual (bs, tmp)) {
      norm = TRUE;
    }
    BSFree (bs);
    bs = tmp;

    SeriousSeqEntryCleanup (sep, NULL, NULL);
    RemoveAllNcbiCleanupUserObjects (sep);
    NormalizeDescriptorOrder (sep);
    tmp = Se2Bs (sep);
    if (! BSEqual (bs, tmp)) {
      ssec = TRUE;
    }
    BSFree (bs);
    bs = tmp;

    BSFree (bs);

    if (ssec) {
      (cfp->rawcounts.ssec)++;
      (cfp->cumcounts.ssec)++;
      if (cfp->logfp != NULL) {
        fprintf (cfp->logfp, "SSEC %s\n", cfp->buf);
        fflush (cfp->logfp);
      }
    } else if (norm) {
      (cfp->rawcounts.norm)++;
      (cfp->cumcounts.norm)++;
      if (cfp->logfp != NULL) {
        fprintf (cfp->logfp, "NORM %s\n", cfp->buf);
        fflush (cfp->logfp);
      }
    } else {
      (cfp->rawcounts.okay)++;
      (cfp->cumcounts.okay)++;
      if (cfp->logfp != NULL) {
        fprintf (cfp->logfp, "OKAY %s\n", cfp->buf);
        fflush (cfp->logfp);
      }
    }

    return;
  }

  MemSet ((Pointer) &cdbefore, 0, sizeof (ChangeData));
  MemSet ((Pointer) &cdafter, 0, sizeof (ChangeData));

  cdbefore.isRefSeq = cfp->isRefSeq;
  cdafter.isRefSeq = cfp->isRefSeq;

  CheckForChanges (sep, &cdbefore);

  bs = Se2Bs (sep);

  NormalizeDescriptorOrder (sep);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    norm = TRUE;
  }
  BSFree (bs);
  bs = tmp;

  VisitFeaturesInSep (sep, NULL, CleanupLocation);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    sloc = TRUE;
  }
  BSFree (bs);
  bs = tmp;

  VisitFeaturesInSep (sep, NULL, CleanupMostRNAs);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    clnr = TRUE;
  }
  BSFree (bs);
  bs = tmp;

  VisitFeaturesInSep (sep, NULL, CleanupRemainingRNAs);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    othr = TRUE;
  }
  BSFree (bs);
  bs = tmp;

  VisitFeaturesInSep (sep, NULL, ModRNAs);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    modr = TRUE;
  }
  BSFree (bs);
  bs = tmp;

  VisitPubdescsInSep (sep, NULL, CleanupPubAuthors);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    auth = TRUE;
  }
  BSFree (bs);
  bs = tmp;

  VisitPubdescsInSep (sep, (Pointer) cfp, CleanupPubBody);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    publ = TRUE;
  }
  BSFree (bs);
  bs = tmp;

  SortSeqEntryQualifiers (sep);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    sort = TRUE;
  }
  BSFree (bs);
  bs = tmp;

  VisitBioseqsInSep (sep, NULL, FixBspMol);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    dnarna = TRUE;
  }
  BSFree (bs);
  bs = tmp;

  BasicSeqEntryCleanup (sep);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    bsec = TRUE;
  }
  BSFree (bs);
  bs = tmp;

  EntryChangeGBSource (sep);
  EntryCheckGBBlock (sep);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    gbbk = TRUE;
  }
  BSFree (bs);
  bs = tmp;

  entityID = ObjMgrGetEntityIDForChoice (sep);
  SeqMgrIndexFeatures (entityID, NULL);
  SeqEntryExplore (sep, NULL, BadProtTitleProc);
  DeleteMarkedObjects (0, OBJ_SEQENTRY, (Pointer) sep);
  SeqMgrIndexFeatures (entityID, NULL);
  InstantiateProteinTitles (entityID, NULL);
  SeqMgrClearFeatureIndexes (entityID, NULL);
  BasicSeqEntryCleanup (sep);
  NormalizeDescriptorOrder (sep);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    titl = TRUE;
  }
  BSFree (bs);
  bs = tmp;

  MoveFeatsFromPartsSet (sep);
  move_cds_ex (sep, TRUE);
  NormalizeDescriptorOrder (sep);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    pack = TRUE;
  }
  BSFree (bs);
  bs = tmp;

  SeqEntryExplore(sep, NULL, ChkSegset);
  NormalizeDescriptorOrder (sep);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    chkseg = TRUE;
  }
  BSFree (bs);
  bs = tmp;

  SeqEntryExplore(sep, (Pointer) &chk_nuc_prot_val, ChkNucProt);
  NormalizeDescriptorOrder (sep);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    chknps = TRUE;
  }
  BSFree (bs);
  bs = tmp;

  SeqEntryPubsAsn4 (sep, cfp->isEmblDdbj);
  NormalizeDescriptorOrder (sep);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    move = TRUE;
  }
  BSFree (bs);
  bs = tmp;

  if (dossec) {
    SeriousSeqEntryCleanup (sep, NULL, NULL);
    RemoveAllNcbiCleanupUserObjects (sep);
    NormalizeDescriptorOrder (sep);
    tmp = Se2Bs (sep);
    if (! BSEqual (bs, tmp)) {
      ssec = TRUE;
    }
    BSFree (bs);
    bs = tmp;
  }

  BSFree (bs);

  CheckForChanges (sep, &cdafter);

  if (ssec || chkseg || chknps) {
    (cfp->rawcounts.ssec)++;
    (cfp->cumcounts.ssec)++;
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "SSEC %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  } else if (move) {
    (cfp->rawcounts.move)++;
    (cfp->cumcounts.move)++;
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "MOVE %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  } else if (pack) {
    (cfp->rawcounts.pack)++;
    (cfp->cumcounts.pack)++;
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "PACK %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  } else if (titl) {
    (cfp->rawcounts.titl)++;
    (cfp->cumcounts.titl)++;
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "TITL %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  } else if (gbbk) {
    (cfp->rawcounts.gbbk)++;
    (cfp->cumcounts.gbbk)++;
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "GBBK %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  } else if (bsec || dnarna) {
    (cfp->rawcounts.bsec)++;
    (cfp->cumcounts.bsec)++;
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "BSEC %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  } else if (sort) {
    (cfp->rawcounts.sort)++;
    (cfp->cumcounts.sort)++;
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "SORT %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  } else if (sloc) {
    (cfp->rawcounts.sloc)++;
    (cfp->cumcounts.sloc)++;
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "SLOC %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  } else if (clnr) {
    (cfp->rawcounts.clnr)++;
    (cfp->cumcounts.clnr)++;
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "CLNR %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  } else if (othr) {
    (cfp->rawcounts.othr)++;
    (cfp->cumcounts.othr)++;
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "OTHR %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  } else if (modr) {
    (cfp->rawcounts.modr)++;
    (cfp->cumcounts.modr)++;
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "MODR %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  } else if (publ) {
    (cfp->rawcounts.publ)++;
    (cfp->cumcounts.publ)++;
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "PUBL %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  } else if (auth) {
    (cfp->rawcounts.auth)++;
    (cfp->cumcounts.auth)++;
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "AUTH %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  } else if (norm) {
    (cfp->rawcounts.norm)++;
    (cfp->cumcounts.norm)++;
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "NORM %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  } else {
    (cfp->rawcounts.okay)++;
    (cfp->cumcounts.okay)++;
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "OKAY %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }

  if (cdbefore.oldgbqual) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "GBQ %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.sgml) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "SGM %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.cdscodon) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "CDN %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.rubisco) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "RUB %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.rbc) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "RBC %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.its) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "ITS %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.rnaother) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "RNA %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.trnanote) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "TRN %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.oldbiomol) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "MOL %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.badOrg) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "ORG %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.rpt_unit_seq) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "RUS %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }

  if (cdbefore.badDbxref) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "BDX %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.refDbxref) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "FDX %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.srcDbxref) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "SDX %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.capDbxref) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "CDX %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.privDbxref) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "PDX %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.oldDbxref) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "ODX %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.multDbxref) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "MDX %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.rareDbxref) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "RDX %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdafter.hasUnpublished && ! cdafter.hasPublished) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "UNP %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }

  if (sort) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "SRT %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (sloc) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "SLC %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (clnr) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "RCN %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (othr) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "RNO %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (modr) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "RMD %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (publ) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "PBC %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (auth) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "ATH %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (pack) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "PKG %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (move) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "MVP %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (titl) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "TTL %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (chkseg) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "SEG %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (chknps) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "NSS%d %s\n", (int) chk_nuc_prot_val, cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (dnarna) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "MRN %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }

  if (cdbefore.protdesc != cdafter.protdesc) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "PRT %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.sfpnote != cdafter.sfpnote) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "COM %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.gbsource != cdafter.gbsource) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "SRC %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  if (cdbefore.cdsconf != cdafter.cdsconf) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "CNF %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }

  if (ncbiusrobj) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "USR %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
}

static CharPtr ffmod [] = {
  "",
  "release",
  "entrez",
  "gbench",
  "dump",
  NULL
};

static void DoGBFFReport (
  SeqEntryPtr sep,
  CleanFlagPtr cfp,
  Int2 batch
)

{
#ifdef OS_UNIX
  AsnIoPtr     aip;
  Char         arguments [128];
  BioseqPtr    bsp;
  Char         ch;
  Char         cmmd [512];
  int          diff;
  FileCache    fc;
  FILE         *fp;
  SeqEntryPtr  fsep;
  Char         line [512];
  FILE         *ofp;
  Char         path1 [PATH_MAX];
  Char         path2 [PATH_MAX];
  Char         path3 [PATH_MAX];
  CharPtr      rep = "reports";
  SeqIdPtr     sip;
  CharPtr      str;

  if (sep == NULL || cfp == NULL) return;

  fsep = FindNthBioseq (sep, 1);
  if (fsep != NULL && fsep->choice == 1) {
    bsp = (BioseqPtr) fsep->data.ptrvalue;
    if (bsp != NULL) {
      for (sip = bsp->id; sip != NULL; sip = sip->next) {
        switch (sip->choice) {
          case SEQID_GENBANK :
            rep = "gbreports";
            break;
          case SEQID_EMBL :
            rep = "ebreports";
            break;
          case SEQID_DDBJ :
            rep = "djreports";
            break;
          case SEQID_OTHER :
            rep = "rfreports";
            break;
          default :
            break;
        }
      }
    }
  }

  if (cfp->logfp != NULL) {
    fprintf (cfp->logfp, "%s\n", cfp->buf);
    fflush (cfp->logfp);
  }

  if (batch == 1) {

    TmpNam (path1);
    TmpNam (path2);

    fp = FileOpen (path1, "w");
    if (fp != NULL) {
      SeqEntryToGnbk (sep, NULL, GENBANK_FMT, cfp->ffmode, NORMAL_STYLE, 0, 0, 0, NULL, fp);
    }
    FileClose (fp);
    SeriousSeqEntryCleanupBulk (sep);
    fp = FileOpen (path2, "w");
    if (fp != NULL) {
      SeqEntryToGnbk (sep, NULL, GENBANK_FMT, cfp->ffmode, NORMAL_STYLE, 0, 0, 0, NULL, fp);
    }
    FileClose (fp);

    sprintf (cmmd, "%s -o %s -n %s -d %s", cfp->ffdiff, path1, path2, rep);
    system (cmmd);

    sprintf (cmmd, "rm %s; rm %s", path1, path2);
    system (cmmd);

  } else if (batch == 2) {

    TmpNam (path1);
    TmpNam (path2);
    TmpNam (path3);

    aip = AsnIoOpen (path3, "w");
    if (aip == NULL) return;

    SeqEntryAsnWrite (sep, aip, NULL);
    AsnIoClose (aip);

    fp = FileOpen (path1, "w");
    if (fp != NULL) {
      SeqEntryToGnbk (sep, NULL, GENBANK_FMT, cfp->ffmode, NORMAL_STYLE, 0, 0, 0, NULL, fp);
    }
    FileClose (fp);

    arguments [0] = '\0';
    sprintf (arguments,
             "-format genbank -mode %s -style normal -view nuc -nocleanup",
             ffmod [(int) cfp->ffmode]);

    sprintf (cmmd, "%s %s -i %s -o %s", cfp->asn2flat, arguments, path3, path2);
    system (cmmd);

    sprintf (cmmd, "diff -h %s %s > %s", path1, path2, path3);
    diff = system (cmmd);

    if (diff > 0) {
      fp = FileOpen (path3, "r");
      ofp = FileOpen (rep, "a");
      if (fp != NULL && ofp != NULL) {
        fprintf (ofp, "\n\n%s\n", cfp->buf);
        fflush (ofp);
        if (FileCacheSetup (&fc, fp)) {
          str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
          while (str != NULL) {
            ch = line [0];
            if (ch == '<' || ch == '>' || ch == '-') {
              fprintf (ofp, "%s\n", line);
            } else if (IS_DIGIT (ch)) {
              fprintf (ofp, "\n%s\n", "===");
            }
            str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
          }
        }
      }
      FileClose (ofp);
      FileClose (fp);
    }

    sprintf (cmmd, "rm %s; rm %s; rm %s", path1, path2, path3);
    system (cmmd);
  }
#endif
}

static void DoValidatorReport (
  SeqEntryPtr sep,
  FILE *logfp,
  CharPtr id,
  CharPtr asnval
)

{
#ifdef OS_UNIX
  AsnIoPtr   aip;
  Char       ch;
  Char       cmmd [512];
  int        diff;
  FileCache  fc;
  FILE       *fp;
  Char       line [512];
  Char       path1 [PATH_MAX];
  Char       path2 [PATH_MAX];
  Char       path3 [PATH_MAX];
  Char       path4 [PATH_MAX];
  Char       path5 [PATH_MAX];
  CharPtr    str;

  if (sep == NULL || logfp == NULL) return;
  if (StringHasNoText (id) || StringHasNoText (asnval)) return;

  TmpNam (path1);
  TmpNam (path2);
  TmpNam (path3);
  TmpNam (path4);
  TmpNam (path5);

  RemoveAllNcbiCleanupUserObjects (sep);

  aip = AsnIoOpen (path3, "w");
  if (aip == NULL) return;

  SeqEntryAsnWrite (sep, aip, NULL);
  AsnIoClose (aip);

  SeriousSeqEntryCleanup (sep, NULL, NULL);
  RemoveAllNcbiCleanupUserObjects (sep);

  aip = AsnIoOpen (path4, "w");
  if (aip == NULL) return;

  SeqEntryAsnWrite (sep, aip, NULL);
  AsnIoClose (aip);

  sprintf (cmmd, "%s -i %s -o stdout -Q 1 -r -l | sort > %s", asnval, path3, path1);
  system (cmmd);

  sprintf (cmmd, "%s -i %s -o stdout -Q 1 -r -l | sort > %s", asnval, path4, path2);
  system (cmmd);

  sprintf (cmmd, "diff -h %s %s > %s", path1, path2, path5);
  diff = system (cmmd);

  if (diff > 0) {
    fp = FileOpen (path5, "r");
    if (fp != NULL) {
      fprintf (logfp, "\n\n%s\n", id);
      fflush (logfp);
      if (FileCacheSetup (&fc, fp)) {
        str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
        while (str != NULL) {
          ch = line [0];
          if (ch == '<' || ch == '>' || ch == '-') {
            fprintf (logfp, "%s\n", line);
          } else if (IS_DIGIT (ch)) {
            fprintf (logfp, "\n%s\n", "===");
          }
          str = FileCacheReadLine (&fc, line, sizeof (line), NULL);
        }
      }
      fflush (logfp);
    }
    FileClose (fp);
  }

  sprintf (cmmd, "rm %s; rm %s; rm %s; rm %s; rm %s", path1, path2, path3, path4, path5);
  system (cmmd);
#endif
}

static void DoModernizeReport (
  SeqEntryPtr sep,
  CleanFlagPtr cfp
)

{
  ByteStorePtr  bs = NULL, tmp = NULL;

  bs = Se2Bs (sep);

  VisitFeaturesInSep (sep, NULL, ModGenes);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "GEN %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  BSFree (bs);
  bs = tmp;

  VisitFeaturesInSep (sep, NULL, ModRNAs);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "NCR %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  BSFree (bs);
  bs = tmp;

  VisitBioSourcesInSep (sep, NULL, ModPCRs);
  tmp = Se2Bs (sep);
  if (! BSEqual (bs, tmp)) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "PCR %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  }
  BSFree (bs);
  bs = tmp;

  BSFree (bs);
}

static CharPtr GetLocString (
  SeqLocPtr slp
)

{
  return SeqLocPrint (slp);
}

static void DoOverlapReport (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  BioseqPtr          bsp;
  Char               buf0 [1000], buf1 [1000], buf2 [1000];
  CleanFlagPtr       cfp;
  SeqMgrFeatContext  fcontext;
  SeqFeatPtr         feat;
  CharPtr            key, locstr1, locstr2;
  SeqIdPtr           sip, siphead;

  if (sfp == NULL || sfp->location == NULL || userdata == NULL) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp->logfp == NULL) return;

  if (sfp->data.choice == SEQFEAT_GENE) return;

  bsp = BioseqFindFromSeqLoc (sfp->location);
  if (bsp == NULL) return;
  if (! ISA_na (bsp->mol)) return;

  siphead = SeqIdSetDup (bsp->id);
  for (sip = siphead; sip != NULL; sip = sip->next) {
    SeqIdStripLocus (sip);
  }
  SeqIdWrite (siphead, buf0, PRINTID_FASTA_LONG, sizeof (buf0));
  SeqIdSetFree (siphead);

  key = FindKeyFromFeatDefType (sfp->idx.subtype, FALSE);
  if (key == NULL) {
    key = "?";
  }

  locstr1 = GetLocString (sfp->location);
  if (locstr1 == NULL) {
    locstr1 = StringSave ("?");
  }

  buf1 [0] = '\0';
  FeatDefLabel (sfp, buf1, sizeof (buf1) - 1, OM_LABEL_BOTH);

  feat = SeqMgrGetOverlappingFeature (sfp->location, FEATDEF_GENE, NULL, 0, NULL, CONTAINED_WITHIN, &fcontext);
  if (feat != NULL) {
    buf2 [0] = '\0';
    FeatDefLabel (feat, buf2, sizeof (buf2) - 1, OM_LABEL_BOTH);
    locstr2 = GetLocString (feat->location);
    if (locstr2 == NULL) {
      locstr2 = StringSave ("?");
    }
    fprintf (cfp->logfp, "%s\t%s\t%s\t%s\t%s\n", buf0, buf1, locstr1, buf2, locstr2);
    MemFree (locstr2);
  }

  if (sfp->data.choice != SEQFEAT_CDREGION) return;

  feat = SeqMgrGetOverlappingFeature (sfp->location, FEATDEF_mRNA, NULL, 0, NULL, CHECK_INTERVALS, &fcontext);
  if (feat != NULL) {
    buf2 [0] = '\0';
    FeatDefLabel (feat, buf2, sizeof (buf2) - 1, OM_LABEL_BOTH);
    locstr2 = GetLocString (feat->location);
    if (locstr2 == NULL) {
      locstr2 = StringSave ("?");
    }
    fprintf (cfp->logfp, "%s\t%s\t%s\t%s\t%s\n", buf0, buf1, locstr1, buf2, locstr2);
    MemFree (locstr2);
  }

  MemFree (locstr1);
}

static CharPtr stopWords [] = {
  "a",
  "about",
  "again",
  "all",
  "almost",
  "also",
  "although",
  "always",
  "among",
  "an",
  "and",
  "another",
  "any",
  "are",
  "as",
  "at",
  "be",
  "because",
  "been",
  "before",
  "being",
  "between",
  "both",
  "but",
  "by",
  "can",
  "could",
  "did",
  "do",
  "does",
  "done",
  "due",
  "during",
  "each",
  "either",
  "enough",
  "especially",
  "etc",
  "for",
  "found",
  "from",
  "further",
  "had",
  "has",
  "have",
  "having",
  "here",
  "how",
  "however",
  "i",
  "if",
  "in",
  "into",
  "is",
  "it",
  "its",
  "itself",
  "just",
  "kg",
  "km",
  "made",
  "mainly",
  "make",
  "may",
  "mg",
  "might",
  "ml",
  "mm",
  "most",
  "mostly",
  "must",
  "nearly",
  "neither",
  "no",
  "nor",
  "obtained",
  "of",
  "often",
  "on",
  "our",
  "overall",
  "perhaps",
  "pmid",
  "quite",
  "rather",
  "really",
  "regarding",
  "seem",
  "seen",
  "several",
  "should",
  "show",
  "showed",
  "shown",
  "shows",
  "significantly",
  "since",
  "so",
  "some",
  "such",
  "than",
  "that",
  "the",
  "their",
  "theirs",
  "them",
  "then",
  "there",
  "therefore",
  "these",
  "they",
  "this",
  "those",
  "through",
  "thus",
  "to",
  "upon",
  "use",
  "used",
  "using",
  "various",
  "very",
  "was",
  "we",
  "were",
  "what",
  "when",
  "which",
  "while",
  "with",
  "within",
  "without",
  "would",
  NULL
};

static Boolean IsStopWord (
  CharPtr str
)

{
  Int2  i;

  if (StringHasNoText (str)) return FALSE;

  for (i = 0; stopWords [i] != NULL; i++) {
    if (StringICmp (str, stopWords [i]) == 0) return TRUE;
  }

  return FALSE;
}

static ValNodePtr GetAuthorMLNameList (
  AuthListPtr alp
)

{
  AuthorPtr    ap;
  Char         buf [128];
  Char         ch;
  Char         chr [4];
  ValNodePtr   head = NULL;
  Char         initials [32];
  ValNodePtr   last = NULL;
  NameStdPtr   nsp;
  PersonIdPtr  pid;
  CharPtr      ptr;
  CharPtr      str;
  ValNodePtr   tmp;
  ValNodePtr   vnp;

  if (alp == NULL) return NULL;

  for (vnp = alp->names; vnp != NULL; vnp = vnp->next) {
    buf [0] = '\0';
    initials [0] = '\0';
    switch (alp->choice) {
      case 1 :
        ap = (AuthorPtr) vnp->data.ptrvalue;
        if (ap == NULL) continue;
        pid = ap->name;
        if (pid == NULL) continue;
        if (pid->choice == 2) {
          nsp = pid->data;
          if (nsp == NULL) continue;
          str = nsp->names [0];
          if (StringHasNoText (str)) continue;
          StringNCpy_0 (buf, str, sizeof (buf));
          StringNCpy_0 (initials, nsp->names [4], sizeof (initials));
        }
        break;
      case 2 :
      case 3 :
        str = (CharPtr) vnp->data.ptrvalue;
        if (StringHasNoText (str)) continue;
        StringNCpy_0 (buf, str, sizeof (buf));
        ptr = StringChr (buf, ',');
        if (ptr == NULL) {
          ptr = StringChr (buf, ' ');
        }
        if (ptr != NULL) {
          *ptr = '\0';
          ptr++;
          StringNCpy_0 (initials, ptr, sizeof (initials));
        }
        break;
      default :
        break;
    }
    if (StringHasNoText (buf)) continue;
    if (StringDoesHaveText (initials)) {
      StringCat (buf, " ");
      chr [1] = '\0';
      ptr = initials;
      ch = *ptr;
      while (ch != '\0') {
        if (ch != ' ' && ch != '.' && ch != ',') {
          chr [0] = ch;
          StringCat (buf, chr);
        }
        ptr++;
        ch = *ptr;
      }
    }
    TrimSpacesAroundString (buf);
    tmp = ValNodeCopyStr (&last, 0, buf);
    if (head == NULL) {
      head = tmp;
    }
    last = tmp;
  }

  return head;
}

static ValNodePtr GetTitleWords (
  CharPtr title
)

{
  Char        ch;
  Boolean     goOn = TRUE;
  ValNodePtr  head = NULL;
  ValNodePtr  last = NULL;
  CharPtr     ptr;
  CharPtr     str;
  CharPtr     tmp;
  ValNodePtr  vnp;

  if (StringHasNoText (title)) return NULL;

  tmp = StringSave (title);
  if (tmp == NULL) return NULL;

  ptr = tmp;
  ch = *ptr;
  if (ch == '\0') {
    goOn = FALSE;
  }
  while (goOn) {
    while (ch != '\0' && (! IS_ALPHANUM (ch))) {
      ptr++;
      ch = *ptr;
    }
    str = ptr;
    while (ch != '\0' && IS_ALPHANUM (ch)) {
      ptr++;
      ch = *ptr;
    }
    if (ch == '\0') {
      goOn = FALSE;
    }
    *ptr = '\0';
    ptr++;
    ch = *ptr;
    TrimSpacesAroundString (str);
    /*
    if (! IsStopWord (str)) {
      vnp = ValNodeCopyStr (&last, 0, str);
      if (head == NULL) {
        head = vnp;
      }
      last = vnp;
    }
    */
    vnp = ValNodeCopyStr (&last, 0, str);
    if (head == NULL) {
      head = vnp;
    }
    last = vnp;
  }

  MemFree (tmp);

  return head;
}

static ValNodePtr DuplicateStringList (
  ValNodePtr list
)

{
  ValNodePtr  head = NULL;
  ValNodePtr  last = NULL;
  CharPtr     str;
  ValNodePtr  tmp;
  ValNodePtr  vnp;

  if (list == NULL) return NULL;

  for (vnp = list; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) continue;
    tmp = ValNodeCopyStr (&last, 0, str);
    if (head == NULL) {
      head = tmp;
    }
    last = tmp;
  }

  return head;
}

typedef enum {
  FULL_INITIALS,
  TWO_INITIALS,
  ONE_INITIAL,
  NO_INITIALS
} InitialsPolicy;

static void TrimInitials (
  CharPtr auth,
  InitialsPolicy initials
)

{
  Char     ch;
  CharPtr  ptr;

  if (StringHasNoText (auth)) return;

  switch (initials) {
    case FULL_INITIALS :
      break;
    case TWO_INITIALS :
      ptr = StringRChr (auth, ' ');
      if (ptr != NULL) {
        ptr++;
        ch = *ptr;
        if (IS_ALPHANUM (ch)) {
          ptr++;
          ch = *ptr;
          if (IS_ALPHANUM (ch)) {
            ptr++;
            *ptr = '\0';
          }
        }
      }
      break;
    case ONE_INITIAL :
      ptr = StringRChr (auth, ' ');
      if (ptr != NULL) {
        ptr++;
        ch = *ptr;
        if (IS_ALPHANUM (ch)) {
          ptr++;
          *ptr = '\0';
        }
      }
      break;
    case NO_INITIALS :
      ptr = StringRChr (auth, ' ');
      if (ptr != NULL) {
        *ptr = '\0';
      }
      break;
    default :
      break;
  }
}

static Int4 DoUnpubBooleanQuery (
  ValNodePtr authors,
  InitialsPolicy initials,
  Boolean firstLastOnly,
  ValNodePtr titlewords,
  Int2 year,
  Boolean expand,
  Uint4Ptr uidp
)

{
  Boolean                 addOpAnd = FALSE;
  Char                    buf [128];
  Int4                    count = 0;
  Entrez2BooleanReplyPtr  e2br;
  Entrez2IdListPtr        e2id;
  Entrez2RequestPtr       e2rp = NULL;
  Entrez2ReplyPtr         e2ry;
  CharPtr                 str;
  ValNodePtr              vnp;

  if (uidp != NULL) {
    *uidp = 0;
  }

  e2rp = EntrezCreateBooleanRequest (TRUE, FALSE, "PubMed", NULL, 0, 0, NULL, 20, 0);
  if (e2rp == NULL) return 0;

  for (vnp = authors; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) continue;
    if (firstLastOnly) {
      if (vnp != authors && vnp->next != NULL) continue;
    }
    StringNCpy_0 (buf, str, sizeof (buf));
    switch (initials) {
      case TWO_INITIALS :
        TrimInitials (buf, TWO_INITIALS);
        break;
      case ONE_INITIAL :
        TrimInitials (buf, ONE_INITIAL);
        break;
      case NO_INITIALS :
        TrimInitials (buf, NO_INITIALS);
        break;
      default :
        break;
    }
    if (addOpAnd) {
      EntrezAddToBooleanRequest (e2rp, NULL, ENTREZ_OP_AND, NULL, NULL, NULL, 0, 0, NULL, NULL, FALSE, FALSE);
    }
    EntrezAddToBooleanRequest (e2rp, NULL, 0, "AUTH", buf, NULL, 0, 0, NULL, NULL, FALSE, FALSE);
    addOpAnd = TRUE;
  }

  for (vnp = titlewords; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) continue;
    if (IsStopWord (str)) continue;
    StringNCpy_0 (buf, str, sizeof (buf));
    if (addOpAnd) {
      EntrezAddToBooleanRequest (e2rp, NULL, ENTREZ_OP_AND, NULL, NULL, NULL, 0, 0, NULL, NULL, FALSE, FALSE);
    }
    EntrezAddToBooleanRequest (e2rp, NULL, 0, "TITL", buf, NULL, 0, 0, NULL, NULL, FALSE, FALSE);
    addOpAnd = TRUE;
  }

  if (year > 0) {
    if (addOpAnd) {
      EntrezAddToBooleanRequest (e2rp, NULL, ENTREZ_OP_AND, NULL, NULL, NULL, 0, 0, NULL, NULL, FALSE, FALSE);
    }
    if (expand) {
      EntrezAddToBooleanRequest (e2rp, NULL, ENTREZ_OP_LEFT_PAREN, NULL, NULL, NULL, 0, 0, NULL, NULL, FALSE, FALSE);
      sprintf (buf, "%d", (int) year - 1);
      EntrezAddToBooleanRequest (e2rp, NULL, 0, "EDAT", buf, NULL, 0, 0, NULL, NULL, FALSE, FALSE);
      EntrezAddToBooleanRequest (e2rp, NULL, ENTREZ_OP_OR, NULL, NULL, NULL, 0, 0, NULL, NULL, FALSE, FALSE);
      sprintf (buf, "%d", (int) year);
      EntrezAddToBooleanRequest (e2rp, NULL, 0, "EDAT", buf, NULL, 0, 0, NULL, NULL, FALSE, FALSE);
      EntrezAddToBooleanRequest (e2rp, NULL, ENTREZ_OP_OR, NULL, NULL, NULL, 0, 0, NULL, NULL, FALSE, FALSE);
      sprintf (buf, "%d", (int) year + 1);
      EntrezAddToBooleanRequest (e2rp, NULL, 0, "EDAT", buf, NULL, 0, 0, NULL, NULL, FALSE, FALSE);
      EntrezAddToBooleanRequest (e2rp, NULL, ENTREZ_OP_RIGHT_PAREN, NULL, NULL, NULL, 0, 0, NULL, NULL, FALSE, FALSE);
    } else {
      sprintf (buf, "%d", (int) year);
      EntrezAddToBooleanRequest (e2rp, NULL, 0, "EDAT", buf, NULL, 0, 0, NULL, NULL, FALSE, FALSE);
    }
    addOpAnd = TRUE;
  }

  e2ry = EntrezSynchronousQuery (e2rp);

  e2rp = Entrez2RequestFree (e2rp);
  if (e2ry == NULL) return 0;
  e2br = EntrezExtractBooleanReply (e2ry);
  if (e2br == NULL) return 0;

  count = e2br->count;

  if (count > 0 && uidp != NULL) {
    e2id = e2br->uids;
    if (e2id != NULL && e2id->num == 1 && e2id->uids != NULL) {
      BSSeek (e2id->uids, 0, SEEK_SET);
      *uidp = Nlm_BSGetUint4 (e2id->uids);
    }
  }

  Entrez2BooleanReplyFree (e2br);

  return count;
}

static CharPtr GetBestJournal (
  ValNodePtr journaltitle
)

{
  CharPtr     str;
  ValNodePtr  vnp;

  if (journaltitle == NULL) return NULL;

  for (vnp = journaltitle; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == Cit_title_iso_jta) {
      str = (CharPtr) vnp->data.ptrvalue;
      if (StringHasNoText (str)) continue;
      return str;
    }
  }

  for (vnp = journaltitle; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == Cit_title_name || vnp->choice == Cit_title_jta) {
      str = (CharPtr) vnp->data.ptrvalue;
      if (StringHasNoText (str)) continue;
      return str;
    }
  }

  return NULL;
}

typedef struct pubref {
  ValNodePtr  authors;
  ValNodePtr  titlewords;
  CharPtr     fulltitle;
  CharPtr     uniquestr;
  CharPtr     journal;
  ImprintPtr  imp;
  Int2        year;
  Uint4       pmid;
} PubRef, PNTR PubRefPtr;

static void PrintPubAuthors (
  CleanFlagPtr cfp,
  PubRefPtr prp
)

{
  CharPtr     prefix = "";
  CharPtr     str;
  ValNodePtr  vnp;

  if (cfp == NULL || cfp->logfp == NULL || prp == NULL) return;

  for (vnp = prp->authors; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) continue;
    fprintf (cfp->logfp, "%s%s", prefix, str);
    prefix = ", ";
  }
}

static void PrintPubTitle (
  CleanFlagPtr cfp,
  PubRefPtr prp
)

{
  CharPtr     prefix = "";
  CharPtr     str;
  ValNodePtr  vnp;

  if (cfp == NULL || cfp->logfp == NULL || prp == NULL) return;

  if (StringDoesHaveText (prp->fulltitle)) {
    fprintf (cfp->logfp, "%s%s", prefix, prp->fulltitle);
  } else {
    for (vnp = prp->titlewords; vnp != NULL; vnp = vnp->next) {
      str = (CharPtr) vnp->data.ptrvalue;
      if (StringHasNoText (str)) continue;
      fprintf (cfp->logfp, "%s%s", prefix, str);
      prefix = " ";
    }
  }
}

static void PrintPubJournal (
  CleanFlagPtr cfp,
  PubRefPtr prp
)

{
  DatePtr     dp = NULL;
  ImprintPtr  imp;
  CharPtr     prefix = "";
  Int2        year;

  if (cfp == NULL || cfp->logfp == NULL || prp == NULL) return;

  if (StringHasNoText (prp->journal) && prp->imp == NULL) {
    fprintf (cfp->logfp, "Unpublished");
    prefix = " ";
    if (prp->year > 0) {
      fprintf (cfp->logfp, "%s[%d]", prefix, (int) prp->year);
      prefix = " ";
    }
    return;
  }

  if (StringDoesHaveText (prp->journal)) {
    fprintf (cfp->logfp, "%s%s", prefix, prp->journal);
    prefix = " ";
  }

  imp = prp->imp;
  if (imp != NULL) {
    dp = imp->date;
    if (dp != NULL && dp->data [0] == 1) {
      year = (Int2) dp->data [1] + 1900;
      fprintf (cfp->logfp, "%s[%d]", prefix, (int) year);
      prefix = " ";
    }
    if (StringDoesHaveText (imp->volume)) {
      fprintf (cfp->logfp, "%s%s", prefix, imp->volume);
      prefix = " ";
    }
    /*
    if (StringDoesHaveText (imp->issue)) {
      fprintf (cfp->logfp, "%s(%s)", prefix, imp->issue);
      prefix = " ";
    }
    */
    if (StringDoesHaveText (imp->pages)) {
      fprintf (cfp->logfp, "%s: %s", prefix, imp->pages);
      prefix = " ";
    }
  }

  if (prp->pmid > 0) {
    fprintf (cfp->logfp, "%s<%ld>", prefix, (long) prp->pmid);
    prefix = " ";
  }
}

/*
static void PrintAuthTitle (
  CleanFlagPtr cfp,
  CharPtr label,
  Boolean multiline,
  PubRefPtr prp
)

{
  DatePtr     dp = NULL;
  ImprintPtr  imp;
  CharPtr     prefix;
  CharPtr     separator = " ";
  CharPtr     str;
  ValNodePtr  vnp;
  Int2        year;

  if (cfp == NULL || cfp->logfp == NULL || prp == NULL) return;

  if (multiline) {
    separator = "\n";
  }

  if (StringDoesHaveText (label)) {
    fprintf (cfp->logfp, "%s", label);
    fprintf (cfp->logfp, separator);
  }

  prefix = "Author:  ";
  for (vnp = prp->authors; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) continue;
    fprintf (cfp->logfp, "%s%s", prefix, str);
    prefix = ", ";
  }
  fprintf (cfp->logfp, separator);
  prefix = "Title:   ";
  if (StringDoesHaveText (prp->fulltitle)) {
    fprintf (cfp->logfp, "%s%s", prefix, prp->fulltitle);
  } else {
    for (vnp = prp->titlewords; vnp != NULL; vnp = vnp->next) {
      str = (CharPtr) vnp->data.ptrvalue;
      if (StringHasNoText (str)) continue;
      fprintf (cfp->logfp, "%s%s", prefix, str);
      prefix = " ";
    }
  }
  fprintf (cfp->logfp, separator);
  year = prp->year;
  if (year > 0) {
    fprintf (cfp->logfp, "Year:    %d", (int) year);
    fprintf (cfp->logfp, separator);
  }

  if (StringDoesHaveText (prp->journal)) {
    fprintf (cfp->logfp, "Journal: %s", prp->journal);
    imp = prp->imp;
    if (imp != NULL) {
      dp = imp->date;
      if (dp != NULL && dp->data [0] == 1) {
        year = (Int2) dp->data [1] + 1900;
        fprintf (cfp->logfp, ". [%d]", (int) year);
      }
      if (StringDoesHaveText (imp->volume)) {
        fprintf (cfp->logfp, ". %s", imp->volume);
      }
      if (StringDoesHaveText (imp->issue)) {
        fprintf (cfp->logfp, " (%s)", imp->issue);
      }
      if (StringDoesHaveText (imp->pages)) {
        fprintf (cfp->logfp, " : %s", imp->pages);
      }
    }
    fprintf (cfp->logfp, separator);
  }

  fprintf (cfp->logfp, "\n");
  fflush (cfp->logfp);
}
*/

typedef enum {
  NO_NAME_MATCH,
  LAST_NAME_MATCH,
  ONE_INIT_MATCH,
  TWO_INIT_MATCH,
  FULL_NAME_MATCH
} AuthComp;

static CharPtr authlabel [] = {
  "AUTH_MISMATCH", "LAST_NAMES", "ONE_INIT", "TWO_INITS", "FULL_NAMES"
};

static AuthComp AuthorCompare (
  CharPtr auth1,
  CharPtr auth2
)

{
  Char  buf1 [128];
  Char  buf2 [128];

  if (StringHasNoText (auth1) || StringHasNoText (auth2)) return NO_NAME_MATCH;

  StringNCpy_0 (buf1, auth1, sizeof (buf1));
  StringNCpy_0 (buf2, auth2, sizeof (buf2));

  if (StringICmp (buf1, buf2) == 0) return FULL_NAME_MATCH;

  TrimInitials (buf1, TWO_INITIALS);
  TrimInitials (buf2, TWO_INITIALS);

  if (StringICmp (buf1, buf2) == 0) return TWO_INIT_MATCH;

  TrimInitials (buf1, ONE_INITIAL);
  TrimInitials (buf2, ONE_INITIAL);

  if (StringICmp (buf1, buf2) == 0) return ONE_INIT_MATCH;

  TrimInitials (buf1, NO_INITIALS);
  TrimInitials (buf2, NO_INITIALS);

  if (StringICmp (buf1, buf2) == 0) return LAST_NAME_MATCH;

  return NO_NAME_MATCH;
}

static AuthComp AuthorsIdentical (
  ValNodePtr oldauthors,
  ValNodePtr newauthors
)

{
  AuthComp    curr, rsult = FULL_NAME_MATCH;
  CharPtr     str1, str2;
  ValNodePtr  vnp1, vnp2;

  if (oldauthors == NULL || newauthors == NULL) return NO_NAME_MATCH;

  for (vnp1 = oldauthors, vnp2 = newauthors;
       vnp1 != NULL && vnp2 != NULL;
       vnp1 = vnp1->next, vnp2 = vnp2->next) {
    str1 = (CharPtr) vnp1->data.ptrvalue;
    str2 = (CharPtr) vnp2->data.ptrvalue;
    curr = AuthorCompare (str1, str2);
    if (curr == NO_NAME_MATCH) return NO_NAME_MATCH;
    if (curr < rsult) {
      rsult = curr;
    }
  }

  if (vnp1 != NULL || vnp2 != NULL) return NO_NAME_MATCH;

  return rsult;
}

static AuthComp AuthorInList (
  CharPtr author,
  ValNodePtr newauthors
)

{
  AuthComp    curr, rsult = FULL_NAME_MATCH;
  CharPtr     str;
  ValNodePtr  vnp;

  if (StringHasNoText (author) || newauthors == NULL) return NO_NAME_MATCH;

  for (vnp = newauthors; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    curr = AuthorCompare (author, str);
    if (curr == NO_NAME_MATCH) continue;
    if (curr < rsult) {
      rsult = curr;
    }
  }

  return rsult;
}

static Boolean WordInList (
  CharPtr word,
  ValNodePtr newtitlewords
)

{
  CharPtr     str;
  ValNodePtr  vnp;

  if (StringHasNoText (word) || newtitlewords == NULL) return NO_NAME_MATCH;

  for (vnp = newtitlewords; vnp != NULL; vnp = vnp->next) {
    str = (CharPtr) vnp->data.ptrvalue;
    if (StringHasNoText (str)) continue;
    if (StringICmp (word, str) == 0) return TRUE;
  }

  return FALSE;
}

static void PrintComparison (
  CleanFlagPtr cfp,
  PubRefPtr oldprp,
  PubRefPtr newprp
)

{
  AuthComp    authcomp, curr, best = FULL_NAME_MATCH;
  Int2        matches, total;
  CharPtr     str, str1, str2;
  Boolean     titlsame;
  ValNodePtr  vnp;

  if (cfp == NULL || cfp->logfp == NULL || oldprp == NULL || newprp == NULL) return;

  authcomp = AuthorsIdentical (oldprp->authors, newprp->authors);
  titlsame = (Boolean) (StringICmp (oldprp->fulltitle, newprp->fulltitle) == 0);

  fprintf (cfp->logfp, "PMID %ld", (long) newprp->pmid);
  fprintf (cfp->logfp, "\t");

  fprintf (cfp->logfp, "%s", cfp->buf);
  fprintf (cfp->logfp, "\t");

  fprintf (cfp->logfp, "REF_COUNT %ld", (long) cfp->unpubcount);
  fprintf (cfp->logfp, "\t");

  fprintf (cfp->logfp, "ORIG_NAMES %ld", (long) ValNodeLen (oldprp->authors));
  fprintf (cfp->logfp, "\t");

  fprintf (cfp->logfp, "ADDL_NAMES %ld", (long) (ValNodeLen (newprp->authors) - ValNodeLen (oldprp->authors)));
  fprintf (cfp->logfp, "\t");

  fprintf (cfp->logfp, "ORIG_WORDS %ld", (long) ValNodeLen (oldprp->titlewords));
  fprintf (cfp->logfp, "\t");

  fprintf (cfp->logfp, "ADDL_WORDS %ld", (long) (ValNodeLen (newprp->titlewords) - ValNodeLen (oldprp->titlewords)));
  fprintf (cfp->logfp, "\t");

  if (StringDoesHaveText (oldprp->uniquestr)) {
    fprintf (cfp->logfp, "UNIQ_CIT %s", oldprp->uniquestr);
  } else {
    fprintf (cfp->logfp, "?");
  }
  fprintf (cfp->logfp, "\t");

  if (authcomp != NO_NAME_MATCH) {
    str = authlabel [(int) authcomp];
    fprintf (cfp->logfp, "AUTHORS_SAME [%s]", str);
    fprintf (cfp->logfp, "\t");
  } else {
    total = ValNodeLen (newprp->authors);

    matches = 0;
    for (vnp = oldprp->authors; vnp != NULL; vnp = vnp->next) {
      str = (CharPtr) vnp->data.ptrvalue;
      if (StringHasNoText (str)) continue;
      curr = AuthorInList (str, newprp->authors);
      if (curr == NO_NAME_MATCH) continue;
      matches++;
      if (curr < best) {
        best = curr;
      }
    }

    str = authlabel [(int) best];
    fprintf (cfp->logfp, "AUTHORS_CHANGED [%s] %d / %d", str, (int) matches, (int) total);
    fprintf (cfp->logfp, "\t");
  }

  if (titlsame) {
    fprintf (cfp->logfp, "TITLE_IDENTICAL");
  } else {
    total = 0;
    for (vnp = newprp->titlewords; vnp != NULL; vnp = vnp->next) {
      str = (CharPtr) vnp->data.ptrvalue;
      if (StringHasNoText (str)) continue;
      if (IsStopWord (str)) continue;
      total++;
    }

    matches = 0;
    for (vnp = oldprp->titlewords; vnp != NULL; vnp = vnp->next) {
      str = (CharPtr) vnp->data.ptrvalue;
      if (StringHasNoText (str)) continue;
      if (IsStopWord (str)) continue;
      if (! WordInList (str, newprp->titlewords)) continue;
      matches++;
    }

    str1 = NULL;
    str2 = NULL;
    vnp = oldprp->titlewords;
    if (vnp != NULL) {
      str1 = (CharPtr) vnp->data.ptrvalue;
    }
    vnp = newprp->titlewords;
    if (vnp != NULL) {
      str2 = (CharPtr) vnp->data.ptrvalue;
    }
    if (str1 != NULL && str2 != NULL && StringCmp (str1, str2) == 0 && total > 0 && matches == total) {
      fprintf (cfp->logfp, "TITLE_SIMILAR %d", (int) matches);
    } else if (total > 0 && matches == total) {
      fprintf (cfp->logfp, "TITLE_ALTERED %d", (int) matches);
    } else {
      fprintf (cfp->logfp, "TITLE_DIFFERS %d / %d", (int) matches, (int) total);
    }
  }
  fprintf (cfp->logfp, "\t");

  fprintf (cfp->logfp, "OLD_AUTH ");
  PrintPubAuthors (cfp, oldprp);
  fprintf (cfp->logfp, "\t");

  fprintf (cfp->logfp, "NEW_AUTH ");
  PrintPubAuthors (cfp, newprp);
  fprintf (cfp->logfp, "\t");

  fprintf (cfp->logfp, "OLD_TITL ");
  PrintPubTitle (cfp, oldprp);
  fprintf (cfp->logfp, "\t");

  fprintf (cfp->logfp, "NEW_TITL ");
  PrintPubTitle (cfp, newprp);
  fprintf (cfp->logfp, "\t");

  fprintf (cfp->logfp, "OLD_JOUR ");
  PrintPubJournal (cfp, oldprp);
  fprintf (cfp->logfp, "\t");

  fprintf (cfp->logfp, "NEW_JOUR ");
  PrintPubJournal (cfp, newprp);

  fprintf (cfp->logfp, "\n");

  /*
  if (identical) {
    PrintAuthTitle (cfp, "EXACT:  ", FALSE, newprp);
  } else {
    PrintAuthTitle (cfp, "BEFORE: ", FALSE, oldprp);
    PrintAuthTitle (cfp, "AFTER:  ", FALSE, newprp);
  }
  */
}

static void StrStripSpaces (
  CharPtr str
)

{
  CharPtr  new_str;

  if (str == NULL) return;

  new_str = str;
  while (*str != '\0') {
    *new_str++ = *str;
    if (*str == ' ' || *str == '\t' || *str == '(') {
      for (str++; *str == ' ' || *str == '\t'; str++) continue;
      if (*str == ')' || *str == ',') {
        new_str--;
      }
    } else {
      str++;
    }
  }
  *new_str = '\0';
}

static void StrStripBrackets (
  CharPtr str
)

{
  size_t  len;

  if (str == NULL) return;

  len = StringLen (str);
  if (len < 2) return;

  if (str [0] == '[') {
    str [0] = ' ';
  }

  if (str [len - 1] == ']') {
    str [len - 1] = ' ';
  }
}

static void PrintPubMedCit (
  CleanFlagPtr cfp,
  Uint4 pmid,
  PubRefPtr oldprp
)

{
  CitArtPtr        cap;
  CitJourPtr       cjp;
  DatePtr          dp;
  ImprintPtr       imp;
  MedlineEntryPtr  mep;
  PubmedEntryPtr   pep;
  PubRef           pr;
  CharPtr          str;
  CharPtr          tmp;
  ValNodePtr       vnp;

  if (cfp == NULL || cfp->logfp == NULL || pmid < 1) return;

  MemSet ((Pointer) &pr, 0, sizeof (PubRef));

  pep = PubMedSynchronousQuery (pmid);
  if (pep == NULL) return;

  mep = (MedlineEntryPtr) pep->medent;
  if (mep != NULL && mep->cit != NULL) {
    cap = mep->cit;
    if (cap != NULL) {
      pr.authors = GetAuthorMLNameList (cap->authors);
      for (vnp = cap->title; vnp != NULL; vnp = vnp->next) {
        if (vnp->choice == Cit_title_name) {
          str = (CharPtr) vnp->data.ptrvalue;
          if (StringHasNoText (str)) continue;
          pr.titlewords = GetTitleWords (str);
          tmp = StringSave (str);
          TrimSpacesAndJunkFromEnds (tmp, TRUE);
          s_RemovePeriodFromEnd (tmp);
          StrStripBrackets (tmp);
          StrStripSpaces (tmp);
          pr.fulltitle = tmp;
        }
      }
      if (cap->from == 1) {
        cjp = (CitJourPtr) cap->fromptr;
        if (cjp != NULL) {
          pr.journal = GetBestJournal (cjp->title);
          imp = cjp->imp;
          pr.imp = imp;
          if (imp != NULL) {
            dp = imp->date;
            if (dp != NULL && dp->data [0] == 1) {
              pr.year = (Int2) dp->data [1] + 1900;
            }
          }
        }
      }
      pr.pmid = pmid;
      if (pr.authors != NULL && pr.titlewords != NULL) {
        /*
        PrintAuthTitle (cfp, "RESULT", TRUE, &pr);
        */
        PrintComparison (cfp, oldprp, &pr);
      }
      ValNodeFreeData (pr.authors);
      ValNodeFreeData (pr.titlewords);
      MemFree (pr.fulltitle);
    }
  }

  pep = PubmedEntryFree (pep);
}

static void TryEntrezQueries (
  CleanFlagPtr cfp,
  PubRefPtr prp
)

{
  Int4   count;
  Uint4  pmid = 0;

  if (cfp == NULL || cfp->logfp == NULL || prp == NULL || prp->authors == NULL) return;

  count = DoUnpubBooleanQuery (prp->authors, ONE_INITIAL, FALSE, prp->titlewords, prp->year, TRUE, &pmid);

  if (count > 1) {
    count = DoUnpubBooleanQuery (prp->authors, TWO_INITIALS, FALSE, prp->titlewords, prp->year, TRUE, &pmid);
  }

  if (count < 1) {

    fprintf (cfp->logfp, "0\t%s\tUNPUB\t", cfp->buf);
    PrintPubAuthors (cfp, prp);
    fprintf (cfp->logfp, "\t");
    PrintPubTitle (cfp, prp);
    fprintf (cfp->logfp, "\t");
    PrintPubJournal (cfp, prp);
    fprintf (cfp->logfp, "\n");

  } else if (count > 1) {

    fprintf (cfp->logfp, "0\t%s\tCOUNT %ld\t", cfp->buf, (long) count);
    PrintPubAuthors (cfp, prp);
    fprintf (cfp->logfp, "\t");
    PrintPubTitle (cfp, prp);
    fprintf (cfp->logfp, "\t");
    PrintPubJournal (cfp, prp);
    fprintf (cfp->logfp, "\n");

  } else {

    PrintPubMedCit (cfp, pmid, prp);
  }

  fflush (cfp->logfp);
}

static void CountUnpubPub (
  PubdescPtr pdp,
  Pointer userdata
)

{
  CleanFlagPtr  cfp;
  CitGenPtr     cgp = NULL;
  Boolean       hasUnpublished = FALSE;
  ValNodePtr    vnp;

  if (pdp == NULL) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp == NULL) return;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == PUB_Gen) {
      cgp = (CitGenPtr) vnp->data.ptrvalue;
      if (cgp != NULL) {
        if (StringICmp (cgp->cit, "Unpublished") == 0) {
          if (StringICmp (cgp->title, "Direct Submission") != 0) {
            hasUnpublished = TRUE;
          }
        }
      }
    } else if (vnp->choice == PUB_Muid || vnp->choice == PUB_PMid) {
      return;
    } else if (vnp->choice == PUB_Article || vnp->choice == PUB_Book || vnp->choice == PUB_Man) {
      return;
    }
  }

  if (! hasUnpublished) return;
  if (cgp == NULL) return;

  (cfp->unpubcount)++;
}

static void ProcessUnpubPub (
  PubdescPtr pdp,
  Pointer userdata
)

{
  Char          buf [521];
  CleanFlagPtr  cfp;
  CitGenPtr     cgp = NULL;
  DatePtr       dp = NULL;
  Boolean       hasUnpublished = FALSE;
  PubRef        pr;
  CharPtr       tmp;
  ValNodePtr    vnp, vnpcgp = NULL;
  Int2          year = 0;

  if (pdp == NULL) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp == NULL) return;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == PUB_Gen) {
      cgp = (CitGenPtr) vnp->data.ptrvalue;
      if (cgp != NULL) {
        if (StringICmp (cgp->cit, "Unpublished") == 0) {
          if (StringICmp (cgp->title, "Direct Submission") != 0) {
            hasUnpublished = TRUE;
            vnpcgp = vnp;
          }
        }
      }
    } else if (vnp->choice == PUB_Muid || vnp->choice == PUB_PMid) {
      return;
    } else if (vnp->choice == PUB_Article || vnp->choice == PUB_Book || vnp->choice == PUB_Man) {
      return;
    }
  }

  if (! hasUnpublished) return;
  if (cgp == NULL) return;

  MemSet ((Pointer) &pr, 0, sizeof (PubRef));

  pr.authors = GetAuthorMLNameList (cgp->authors);
  pr.titlewords = GetTitleWords (cgp->title);
  if (vnpcgp != NULL) {
    if (PubLabelUnique (vnpcgp, buf, sizeof (buf) - 1, OM_LABEL_CONTENT, TRUE) > 0) {
      pr.uniquestr = StringSaveNoNull (buf);
    }
  }

  tmp = StringSave (cgp->title);
  TrimSpacesAndJunkFromEnds (tmp, TRUE);
  s_RemovePeriodFromEnd (tmp);
  StrStripBrackets (tmp);
  StrStripSpaces (tmp);
  pr.fulltitle = tmp;

  pr.journal = GetBestJournal (cgp->journal);
  pr.imp = NULL;

  dp = cgp->date;
  if (dp != NULL && dp->data [0] == 1) {
    year = (Int2) dp->data [1] + 1900;
  }
  if (year == 0) {
    year = cfp->year;
  }
  pr.year = year;
  pr.pmid = 0;

  if (pr.authors != NULL && pr.titlewords != NULL) {
    TryEntrezQueries (cfp, &pr);
  }

  ValNodeFreeData (pr.authors);
  ValNodeFreeData (pr.titlewords);
  MemFree (pr.fulltitle);
  MemFree (pr.uniquestr);
}

static void DoUnpublishedReport (
  SeqEntryPtr sep,
  CleanFlagPtr cfp
)

{
  if (sep == NULL || cfp == NULL) return;

  cfp->unpubcount = 0;
  VisitPubdescsInSep (sep, (Pointer) cfp, CountUnpubPub);
  VisitPubdescsInSep (sep, (Pointer) cfp, ProcessUnpubPub);
}

static void CountPublishedPub (
  PubdescPtr pdp,
  Pointer userdata
)

{
  CleanFlagPtr  cfp;
  CitArtPtr     cap = NULL;
  ValNodePtr    vnp;

  if (pdp == NULL) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp == NULL) return;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == PUB_Article) {
      cap = (CitArtPtr) vnp->data.ptrvalue;
    } else if (vnp->choice == PUB_PMid) {
      return;
    }
  }

  if (cap == NULL) return;

  (cfp->unpubcount)++;
}

static void ProcessPublishedPub (
  PubdescPtr pdp,
  Pointer userdata
)

{
  CleanFlagPtr  cfp;
  CitArtPtr     cap = NULL;
  CitJourPtr    cjp;
  DatePtr       dp = NULL;
  ImprintPtr    imp;
  PubRef        pr;
  CharPtr       str;
  CharPtr       tmp;
  ValNodePtr    vnp;
  Int2          year = 0;

  if (pdp == NULL) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp == NULL) return;

  for (vnp = pdp->pub; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == PUB_Article) {
      cap = (CitArtPtr) vnp->data.ptrvalue;
    } else if (vnp->choice == PUB_PMid) {
      return;
    }
  }

  if (cap == NULL) return;

  MemSet ((Pointer) &pr, 0, sizeof (PubRef));

  pr.authors = GetAuthorMLNameList (cap->authors);
  for (vnp = cap->title; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == Cit_title_name) {
      str = (CharPtr) vnp->data.ptrvalue;
      if (StringHasNoText (str)) continue;
      pr.titlewords = GetTitleWords (str);
      tmp = StringSave (str);
      TrimSpacesAndJunkFromEnds (tmp, TRUE);
      s_RemovePeriodFromEnd (tmp);
      StrStripBrackets (tmp);
      StrStripSpaces (tmp);
      pr.fulltitle = tmp;
    }
  }

  if (cap->from == 1) {
    cjp = (CitJourPtr) cap->fromptr;
    if (cjp != NULL) {
      pr.journal = GetBestJournal (cjp->title);
      imp = cjp->imp;
      pr.imp = imp;
      if (imp != NULL) {
        dp = imp->date;
        if (dp != NULL && dp->data [0] == 1) {
          year = (Int2) dp->data [1] + 1900;
        }
      }
    }
  }

  if (year == 0) {
    year = cfp->year;
  }
  pr.year = year;
  pr.pmid = 0;

  if (pr.authors != NULL && pr.titlewords != NULL) {
    TryEntrezQueries (cfp, &pr);
  }

  ValNodeFreeData (pr.authors);
  ValNodeFreeData (pr.titlewords);
  MemFree (pr.fulltitle);
}

static void DoPublishedReport (
  SeqEntryPtr sep,
  CleanFlagPtr cfp
)

{
  if (sep == NULL || cfp == NULL) return;

  cfp->unpubcount = 0;
  VisitPubdescsInSep (sep, (Pointer) cfp, CountPublishedPub);
  VisitPubdescsInSep (sep, (Pointer) cfp, ProcessPublishedPub);
}

static void RemoveFeatureCitations (
  SeqFeatPtr sfp,
  Pointer userdata
)

{
  if (sfp == NULL || sfp->cit == NULL) return;

  sfp->cit = PubSetFree (sfp->cit);
}

#ifdef OS_UNIX
static SeqEntryPtr CppBasicCleanup (
  SeqEntryPtr sep,
  CleanFlagPtr cfp
)

{
  AsnIoPtr      aip, aop;
  ByteStorePtr  bs1, bs2;
  Char          cmmd [512];
  SeqEntryPtr   csep, nsep;
  Char          path1 [PATH_MAX];
  Char          path2 [PATH_MAX];
  Char          path3 [PATH_MAX];

  if (sep == NULL || cfp == NULL) return NULL;

  VisitFeaturesInSep (sep, NULL, RemoveFeatureCitations);

  TmpNam (path1);
  TmpNam (path2);
  TmpNam (path3);

  aop = AsnIoOpen (path1, "w");
  SeqEntryAsnWrite (sep, aop, NULL);
  AsnIoClose (aop);

  sprintf (cmmd, "%s -i %s | cleanasn -a e -o %s",
           "~/ncbi_cxx/compilers/xCode/build/bin/Debug/test_basic_cleanup",
           path1, path2);
  system (cmmd);

  sprintf (cmmd, "cleanasn -i %s -o %s -K b",
           path1, path3);
  system (cmmd);

  aip = AsnIoOpen (path3, "r");
  csep = SeqEntryAsnRead (aip, NULL);
  AsnIoClose (aip);

  bs1 = Se2Bs (csep);

  aip = AsnIoOpen (path2, "r");
  nsep = SeqEntryAsnRead (aip, NULL);
  AsnIoClose (aip);

  bs2 = Se2Bs (nsep);

  if (nsep == NULL) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "EMPTY %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
  } else if (! BSEqual (bs1, bs2)) {
    if (cfp->logfp != NULL) {
      fprintf (cfp->logfp, "BSEC DIFF %s\n", cfp->buf);
      fflush (cfp->logfp);
    }
    if (cfp->gi > 0) {
      sprintf (cmmd, "echo '' >> ~/Desktop/diffclean.txt");
      system (cmmd);
      sprintf (cmmd, "echo '' >> ~/Desktop/diffclean.txt");
      system (cmmd);
      sprintf (cmmd, "echo '********** gi|%ld **********' >> ~/Desktop/diffclean.txt", (long) cfp->gi);
      system (cmmd);
      sprintf (cmmd, "echo '' >> ~/Desktop/diffclean.txt");
      system (cmmd);
      sprintf (cmmd, "diff %s %s >> ~/Desktop/diffclean.txt", path3, path2);
      system (cmmd);
    }
  }

  BSFree (bs1);
  BSFree (bs2);

  SeqEntryFree (csep);

  sprintf (cmmd, "rm %s; rm %s; rm %s", path1, path2, path3);
  system (cmmd);

  return nsep;
}
#endif

/* now only strips serials for local, general, refseq, and 2+6 genbank ids */
static void CheckForSwissProtIDX (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr     bsp;
  SeqIdPtr      sip;
  BoolPtr       stripSerial;
  TextSeqIdPtr  tsip;

  if (sep == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    if (bsp == NULL) return;
    stripSerial = (BoolPtr) mydata;
    if (stripSerial == NULL) return;
    for (sip = bsp->id; sip != NULL; sip = sip->next) {
      switch (sip->choice) {
        case SEQID_GIBBSQ :
        case SEQID_GIBBMT :
          *stripSerial = FALSE;
          break;
        case SEQID_EMBL :
        case SEQID_PIR :
        case SEQID_SWISSPROT :
        case SEQID_PATENT :
        case SEQID_DDBJ :
        case SEQID_PRF :
        case SEQID_PDB :
        case SEQID_TPE:
        case SEQID_TPD:
        case SEQID_GPIPE:
          *stripSerial = FALSE;
          break;
        case SEQID_GENBANK :
        case SEQID_TPG:
          tsip = (TextSeqIdPtr) sip->data.ptrvalue;
          if (tsip != NULL) {
            if (StringLen (tsip->accession) == 6) {
              *stripSerial = FALSE;
            }
          }
          break;
        case SEQID_NOT_SET :
        case SEQID_LOCAL :
        case SEQID_OTHER :
        case SEQID_GENERAL :
          break;
        default :
          break;
      }
    }
  }
}

static time_t DoCleanup (
  SeqEntryPtr sep,
  Uint2 entityID,
  CleanFlagPtr cfp,
  AsnIoPtr aop,
  AsnTypePtr atp,
  SeqSubmitPtr ssp
)

{
  Boolean              all_digits = TRUE, isDdbj = FALSE, isEmbl = FALSE,
                       isGenBank = FALSE, isNcbi = FALSE, isRefSeq = FALSE,
                       stripSerial = TRUE;
  BioseqPtr            bsp;
  Char                 ch;
  DatePtr              dp;
  SeqEntryPtr          fsep, nsep = NULL;
  Int4                 nucs, prts, pmid = 0;
  CharPtr              ptr;
  SumDataPtr           sdp;
  SeqIdPtr             sip, siphead;
  time_t               starttime, stoptime;
  long int             val;
  SeqDescrPtr          vnp;
  PopSetRetroStatData  stat;

  if (sep == NULL || cfp == NULL) return 0;

  MemSet ((Pointer) &stat, 0, sizeof (PopSetRetroStatData));

  AssignIDsInEntityEx (entityID, 0, NULL, NULL);

  starttime = GetSecs ();

  StringCpy (cfp->buf, "");
  cfp->gi = 0;
  cfp->year = 0;
  cfp->isRefSeq = FALSE;
  cfp->isEmblDdbj = FALSE;

  fsep = FindNthBioseq (sep, 1);
  if (fsep != NULL && fsep->choice == 1) {
    bsp = (BioseqPtr) fsep->data.ptrvalue;
    if (bsp != NULL) {
      siphead = SeqIdSetDup (bsp->id);
      for (sip = siphead; sip != NULL; sip = sip->next) {
        SeqIdStripLocus (sip);
        if (sip->choice == SEQID_GI) {
          cfp->gi = (Int4) sip->data.intvalue;
        } else if (sip->choice == SEQID_GENBANK || sip->choice == SEQID_TPG) {
          isGenBank = TRUE;
          isNcbi = TRUE;
        } else if (sip->choice == SEQID_EMBL || sip->choice == SEQID_TPE) {
          isEmbl = TRUE;
          cfp->isEmblDdbj = TRUE;
        } else if (sip->choice == SEQID_DDBJ || sip->choice == SEQID_TPD) {
          isDdbj = TRUE;
          cfp->isEmblDdbj = TRUE;
        } else if (sip->choice == SEQID_OTHER) {
          isRefSeq = TRUE;
          isNcbi = TRUE;
          cfp->isRefSeq = TRUE;
        }
      }
      SeqIdWrite (siphead, cfp->buf, PRINTID_FASTA_LONG, sizeof (cfp->buf));
      SeqIdSetFree (siphead);
    }
    vnp = GetNextDescriptorUnindexed (bsp, Seq_descr_update_date, NULL);
    if (vnp == NULL) {
      vnp = GetNextDescriptorUnindexed (bsp, Seq_descr_create_date, NULL);
    }
    if (vnp != NULL) {
      dp = (DatePtr) vnp->data.ptrvalue;
      if (dp != NULL && dp->data [0] == 1) {
        cfp->year = (Int2) dp->data [1] + 1900;
      }
    }
  }

  SeqEntryExplore (sep, (Pointer) &stripSerial, CheckForSwissProtIDX);
  cfp->stripSerial = stripSerial;

  if (StringDoesHaveText (cfp->sourcedb)) {
    if (StringChr (cfp->sourcedb, 'g') != NULL) {
      if (! isGenBank) return 0;
    }
    if (StringChr (cfp->sourcedb, 'e') != NULL) {
      if (! isEmbl) return 0;
    }
    if (StringChr (cfp->sourcedb, 'd') != NULL) {
      if (! isDdbj) return 0;
    }
    if (StringChr (cfp->sourcedb, 'r') != NULL) {
      if (! isRefSeq) return 0;
    }
    if (StringChr (cfp->sourcedb, 'n') != NULL) {
      if (! isNcbi) return 0;
    }
    if (StringChr (cfp->sourcedb, 'x') != NULL) {
      if (isEmbl || isDdbj) return 0;
    }
  }

  nucs = VisitSequencesInSep (sep, NULL, VISIT_NUCS, NULL);
  prts = VisitSequencesInSep (sep, NULL, VISIT_PROTS, NULL);
  cfp->rawcounts.nucs += nucs;
  cfp->rawcounts.prts += prts;
  (cfp->rawcounts.recs)++;
  cfp->cumcounts.nucs += nucs;
  cfp->cumcounts.prts += prts;
  (cfp->cumcounts.recs)++;

  sdp = NULL;
  if (isGenBank) {
    sdp = &(cfp->dbsums.genbank);
  } else if (isEmbl) {
    sdp = &(cfp->dbsums.embl);
  } else if (isDdbj) {
    sdp = &(cfp->dbsums.ddbj);
  } else if (isRefSeq) {
    sdp = &(cfp->dbsums.refseq);
  } else {
    sdp = &(cfp->dbsums.other);
  }
  if (sdp != NULL) {
    sdp->nucs += nucs;
    sdp->prts += prts;
    (sdp->recs)++;
  }

  if (StringChr (cfp->report, 'c') != NULL) {
    return 0;
  }
  if (StringChr (cfp->report, 'r') != NULL) {
    DoASNReport (sep, cfp, FALSE, FALSE, FALSE);
    stoptime = GetSecs ();
    return stoptime - starttime;
  }
  if (StringChr (cfp->report, 's') != NULL) {
    DoASNReport (sep, cfp, TRUE, FALSE, FALSE);
    stoptime = GetSecs ();
    return stoptime - starttime;
  }
  if (StringChr (cfp->report, 'n') != NULL) {
    DoASNReport (sep, cfp, TRUE, TRUE, FALSE);
    stoptime = GetSecs ();
    return stoptime - starttime;
  }
  if (StringChr (cfp->report, 'e') != NULL) {
    DoASNReport (sep, cfp, FALSE, FALSE, TRUE);
    stoptime = GetSecs ();
    return stoptime - starttime;
  }
  if (StringChr (cfp->report, 'd') != NULL) {
    DoAsnDiffReport (sep, cfp);
    stoptime = GetSecs ();
    return stoptime - starttime;
  }
  if (StringChr (cfp->report, 'g') != NULL) {
    DoGBFFReport (sep, cfp, 1);
    stoptime = GetSecs ();
    return stoptime - starttime;
  }
  if (StringChr (cfp->report, 'f') != NULL) {
    DoGBFFReport (sep, cfp, 2);
    stoptime = GetSecs ();
    return stoptime - starttime;
  }
  if (StringChr (cfp->report, 'v') != NULL) {
    DoValidatorReport (sep, cfp->logfp, cfp->buf, cfp->asnval);
    stoptime = GetSecs ();
    return stoptime - starttime;
  }
  if (StringChr (cfp->report, 'm') != NULL) {
    DoModernizeReport (sep, cfp);
    stoptime = GetSecs ();
    return stoptime - starttime;
  }
  if (StringChr (cfp->report, 'o') != NULL) {
    SeqMgrIndexFeatures (entityID, 0);
    VisitFeaturesInSep (sep, (Pointer) cfp, DoOverlapReport);
    stoptime = GetSecs ();
    return stoptime - starttime;
  }
  if (StringChr (cfp->report, 'u') != NULL) {
    DoUnpublishedReport (sep, cfp);
    stoptime = GetSecs ();
    return stoptime - starttime;
  }
  if (StringChr (cfp->report, 'p') != NULL) {
    DoPublishedReport (sep, cfp);
    stoptime = GetSecs ();
    return stoptime - starttime;
  }
  if (StringChr (cfp->report, 'x') != NULL) {
    DoASNScan (sep, cfp);
    stoptime = GetSecs ();
    return stoptime - starttime;
  }

  if (StringDoesHaveText (cfp->report)) return 0;

  if (cfp->logfp != NULL) {
    fprintf (cfp->logfp, "%s\n", cfp->buf);
    fflush (cfp->logfp);
  }

  if (StringChr (cfp->clean, 'b') != NULL) {
    BasicSeqEntryCleanup (sep);
  }
#ifdef OS_UNIX
  if (StringChr (cfp->clean, 'p') != NULL) {
    nsep = CppBasicCleanup (sep, cfp);
  }
#endif
  if (StringChr (cfp->clean, 's') != NULL) {
    SeriousSeqEntryCleanup (sep, NULL, NULL);
  }
  if (StringChr (cfp->clean, 'g') != NULL) {
    GpipeSeqEntryCleanup (sep);
  }
  if (StringChr (cfp->clean, 'n') != NULL) {
    NormalizeDescriptorOrder (sep);
  }
  if (StringChr (cfp->clean, 'u') != NULL) {
    RemoveAllNcbiCleanupUserObjects (sep);
  }

  if (StringChr (cfp->modernize, 'g') != NULL) {
    VisitFeaturesInSep (sep, NULL, ModGenes);
  }
  if (StringChr (cfp->modernize, 'r') != NULL) {
    VisitFeaturesInSep (sep, NULL, ModRNAs);
  }
  if (StringChr (cfp->modernize, 'p') != NULL) {
    VisitBioSourcesInSep (sep, NULL, ModPCRs);
  }

  if (cfp->taxon) {
    Taxon3ReplaceOrgInSeqEntry (sep, FALSE);
  }

  if (StringChr (cfp->pub, 'u') != NULL) {
    VisitPubdescsInSep (sep, NULL, LookupPubdesc);
  }
  if (cfp->pub != NULL) {
    ptr = cfp->pub;
    ch = *ptr;
    while (ch != '\0') {
      if (! IS_DIGIT (ch)) {
        all_digits = FALSE;
      }
      ptr++;
      ch = *ptr;
    }
    if (all_digits && sscanf (cfp->pub, "%ld", &val) == 1 && val != 0) {
      pmid = (Int4) val;
      VisitPubdescsInSep (sep, (Pointer) &pmid, ReplaceUnpub);
    }
  }

  if (StringChr (cfp->feat, 'u') != NULL) {
    VisitFeaturesInSep (sep, NULL, RemoveFeatUser);
  }
  if (StringChr (cfp->feat, 'd') != NULL) {
    VisitFeaturesInSep (sep, NULL, RemoveFeatDbxref);
  }
  if (StringChr (cfp->feat, 'r') != NULL) {
    SeqMgrIndexFeatures (entityID, 0);
    VisitFeaturesInSep (sep, NULL, RemoveUnnecGeneXref);
  }

  if (StringChr (cfp->link, 'o') != NULL) {
    SeqMgrIndexFeatures (entityID, 0);
    LinkCDSmRNAbyOverlap (sep);
  }
  if (StringChr (cfp->link, 'p') != NULL) {
    SeqMgrIndexFeatures (entityID, 0);
    LinkCDSmRNAbyProduct (sep);
  }
  if (StringChr (cfp->link, 'r') != NULL) {
    SeqMgrIndexFeatures (entityID, 0);
    ReassignFeatureIDs (sep);
  }
  if (StringChr (cfp->link, 'c') != NULL) {
    ClearFeatureIDs (sep);
  }

  if (StringChr (cfp->desc, 't') != NULL) {
    VisitDescriptorsInSep (sep, NULL, MarkTitles);
    DeleteMarkedObjects (entityID, 0, NULL);
  }
  if (StringChr (cfp->desc, 'n') != NULL) {
    RemoveNucProtSetTitles (sep);
    DeleteMarkedObjects (entityID, 0, NULL);
  }
  if (StringChr (cfp->desc, 'e') != NULL) {
    RemovePopsetTitles (sep);
    DeleteMarkedObjects (entityID, 0, NULL);
  }

  if (StringChr (cfp->mods, 'd') != NULL) {
    SeqMgrIndexFeatures (entityID, 0);
    DoAutoDef (sep, entityID);
  }
  if (StringChr (cfp->mods, 'e') != NULL) {
    SeqMgrIndexFeatures (entityID, 0);
    PopSetAutoDefRetro (sep, &stat);
    LogPopSetTitleResults (cfp->logfp, &stat, sep);
  }

  if (cfp->action_list != NULL) {
    ApplyMacroToSeqEntry (sep, cfp->action_list, NULL, NULL);
  }

  /* normalize order again at end */
  if (StringChr (cfp->clean, 'n') != NULL) {
    NormalizeDescriptorOrder (sep);
  }

  stoptime = GetSecs ();

  if (aop != NULL) {
    if (ssp != NULL) {
      SeqSubmitAsnWrite (ssp, aop, atp);
    } else if (nsep != NULL) {
      SeqEntryAsnWrite (nsep, aop, atp);
      SeqEntryFree (nsep);
    } else {
      SeqEntryAsnWrite (sep, aop, atp);
    }
  }

  return stoptime - starttime;
}

static void CleanupSingleRecord (
  CharPtr filename,
  CleanFlagPtr cfp
)

{
  AsnIoPtr      aip, aop = NULL;
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Pointer       dataptr = NULL;
  Uint2         datatype, entityID = 0;
  FILE          *fp;
  Char          path [PATH_MAX];
  CharPtr       ptr;
  SeqEntryPtr   sep;
  SeqSubmitPtr  ssp = NULL;

  if (cfp == NULL) return;

  if (StringHasNoText (filename)) return;

  if (cfp->type == 1) {
    fp = FileOpen (filename, "r");
    if (fp == NULL) {
      Message (MSG_POSTERR, "Failed to open '%s'", filename);
      return;
    }

    dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, NULL, FALSE, FALSE, FALSE, FALSE);

    FileClose (fp);

    entityID = ObjMgrRegister (datatype, dataptr);

  } else if (cfp->type >= 2 && cfp->type <= 5) {
    aip = AsnIoOpen (filename, cfp->binary? "rb" : "r");
    if (aip == NULL) {
      Message (MSG_POSTERR, "AsnIoOpen failed for input file '%s'", filename);
      return;
    }

    SeqMgrHoldIndexing (TRUE);
    switch (cfp->type) {
      case 2 :
        dataptr = (Pointer) SeqEntryAsnRead (aip, NULL);
        datatype = OBJ_SEQENTRY;
        break;
      case 3 :
        dataptr = (Pointer) BioseqAsnRead (aip, NULL);
        datatype = OBJ_BIOSEQ;
        break;
      case 4 :
        dataptr = (Pointer) BioseqSetAsnRead (aip, NULL);
        datatype = OBJ_BIOSEQSET;
        break;
      case 5 :
        dataptr = (Pointer) SeqSubmitAsnRead (aip, NULL);
        ssp = (SeqSubmitPtr) dataptr;
        datatype = OBJ_SEQSUB;
        break;
      default :
        break;
    }
    SeqMgrHoldIndexing (FALSE);

    AsnIoClose (aip);

    entityID = ObjMgrRegister (datatype, dataptr);

  } else {
    Message (MSG_POSTERR, "Input format type '%d' unrecognized", (int) cfp->type);
    return;
  }

  if (entityID < 1 || dataptr == NULL) {
    Message (MSG_POSTERR, "Data read failed for input file '%s'", filename);
    return;
  }

  if (datatype == OBJ_SEQSUB || datatype == OBJ_SEQENTRY ||
        datatype == OBJ_BIOSEQ || datatype == OBJ_BIOSEQSET) {

    sep = GetTopSeqEntryForEntityID (entityID);

    if (sep == NULL) {
      sep = SeqEntryNew ();
      if (sep != NULL) {
        if (datatype == OBJ_BIOSEQ) {
          bsp = (BioseqPtr) dataptr;
          sep->choice = 1;
          sep->data.ptrvalue = bsp;
          SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) bsp, sep);
        } else if (datatype == OBJ_BIOSEQSET) {
          bssp = (BioseqSetPtr) dataptr;
          sep->choice = 2;
          sep->data.ptrvalue = bssp;
          SeqMgrSeqEntry (SM_BIOSEQSET, (Pointer) bssp, sep);
        } else {
          sep = SeqEntryFree (sep);
        }
      }
      sep = GetTopSeqEntryForEntityID (entityID);
    }

    if (sep != NULL) {

      path [0] = '\0';
      if (StringDoesHaveText (cfp->outfile)) {

        StringNCpy_0 (path, cfp->outfile, sizeof (path));

      } else if (StringDoesHaveText (cfp->results)) {

        ptr = StringRChr (filename, DIRDELIMCHR);
        if (ptr != NULL) {
          StringNCpy_0 (path, cfp->results, sizeof (path));
          ptr++;
          FileBuildPath (path, NULL, ptr);
        }
      }

      sep = GetTopSeqEntryForEntityID (entityID);
      if (sep != NULL) {

        if (StringHasNoText (cfp->report) && StringDoesHaveText (path)) {
          aop = AsnIoOpen (path, "w");
        }

        DoCleanup (sep, entityID, cfp, aop, NULL, ssp);

        if (aop != NULL) {
          AsnIoFlush (aop);
          AsnIoClose (aop);
        }
      }

      ObjMgrFreeByEntityID (entityID);
    }

  } else {

    Message (MSG_POSTERR, "Datatype %d not recognized", (int) datatype);
  }
}

static void CleanupMultipleRecord (
  CharPtr filename,
  CleanFlagPtr cfp
)

{
  AsnIoPtr     aip, aop = NULL;
  AsnTypePtr   atp;
  DataVal      av;
  Char         ch;
  Uint2        entityID;
  FILE         *fp;
  size_t       len;
  Char         longest [64];
  Int4         numrecords;
  Char         path [PATH_MAX];
  CharPtr      ptr;
  SeqEntryPtr  sep;
  time_t       timediff, worsttime;
#ifdef OS_UNIX
  Char         cmmd [512];
  CharPtr      gzcatprog;
  int          ret;
  Boolean      usedPopen = FALSE;
#endif

  if (cfp == NULL) return;

  if (StringHasNoText (filename)) return;

  path [0] = '\0';
  if (StringDoesHaveText (cfp->outfile)) {

    StringNCpy_0 (path, cfp->outfile, sizeof (path));

  } else if (StringDoesHaveText (cfp->results)) {

    ptr = StringRChr (filename, DIRDELIMCHR);
    if (ptr != NULL) {
      StringNCpy_0 (path, cfp->results, sizeof (path));
      ptr++;
      if (cfp->compressed) {
        len = StringLen (ptr);
        if (len > 4 && StringCmp (ptr + len - 3, ".gz") == 0) {
          ptr [len - 3] = '\0';
        }
      }
      FileBuildPath (path, NULL, ptr);
    }
  }
  if (StringHasNoText (cfp->report) && StringHasNoText (path)) return;

#ifndef OS_UNIX
  if (cfp->compressed) {
    Message (MSG_POSTERR, "Can only decompress on-the-fly on UNIX machines");
    return;
  }
#endif

#ifdef OS_UNIX
  if (cfp->compressed) {
    gzcatprog = getenv ("NCBI_UNCOMPRESS_BINARY");
    if (gzcatprog != NULL) {
      sprintf (cmmd, "%s %s", gzcatprog, filename);
    } else {
      ret = system ("gzcat -h >/dev/null 2>&1");
      if (ret == 0) {
        sprintf (cmmd, "gzcat %s", filename);
      } else if (ret == -1) {
        Message (MSG_POSTERR, "Unable to fork or exec gzcat in ScanBioseqSetRelease");
        return;
      } else {
        ret = system ("zcat -h >/dev/null 2>&1");
        if (ret == 0) {
          sprintf (cmmd, "zcat %s", filename);
        } else if (ret == -1) {
          Message (MSG_POSTERR, "Unable to fork or exec zcat in ScanBioseqSetRelease");
          return;
        } else {
          Message (MSG_POSTERR, "Unable to find zcat or gzcat in ScanBioseqSetRelease - please edit your PATH environment variable");
          return;
        }
      }
    }
    fp = popen (cmmd, /* cfp->binary? "rb" : */ "r");
    usedPopen = TRUE;
  } else {
    fp = FileOpen (filename, cfp->binary? "rb" : "r");
  }
#else
  fp = FileOpen (filename, cfp->binary? "rb" : "r");
#endif
  if (fp == NULL) {
    Message (MSG_POSTERR, "FileOpen failed for input file '%s'", filename);
    return;
  }

  aip = AsnIoNew (cfp->binary? ASNIO_BIN_IN : ASNIO_TEXT_IN, fp, NULL, NULL, NULL);
  if (aip == NULL) {
    Message (MSG_ERROR, "AsnIoNew failed for input file '%s'", filename);
    return;
  }

  if (cfp->logfp != NULL) {
    if (StringChr (cfp->report, 'c') == NULL) {
      fprintf (cfp->logfp, "%s\n\n", filename);
      fflush (cfp->logfp);
    }
  }

  longest [0] = '\0';
  worsttime = 0;
  numrecords = 0;

  if (StringHasNoText (cfp->report)) {
    aop = AsnIoOpen (path, cfp->binary? "wb" : "w");
    if (aop != NULL) {
      AsnOpenStruct (aop, cfp->bssp_atp, (Pointer) &(cfp->bss));
      av.intvalue = 7;
      AsnWrite (aop, cfp->atp_bsc, &av);
      AsnOpenStruct (aop, cfp->atp_bsss, (Pointer) &(cfp->bss.seq_set));
    }
  }

  atp = cfp->atp_bss;

  while ((atp = AsnReadId (aip, cfp->amp, atp)) != NULL) {
    if (atp == cfp->atp_se) {

      SeqMgrHoldIndexing (TRUE);
      sep = SeqEntryAsnRead (aip, atp);
      SeqMgrHoldIndexing (FALSE);

      if (sep != NULL) {

        entityID = ObjMgrGetEntityIDForChoice (sep);

        timediff = DoCleanup (sep, entityID, cfp, aop, cfp->atp_se, NULL);

        if (timediff > worsttime) {
          worsttime = timediff;
          StringCpy (longest, cfp->buf);
          ptr = longest;
          ch = *ptr;
          while (ch != '\0') {
            if (ch == '|') {
              *ptr = ' ';
            }
            ptr++;
            ch = *ptr;
          }
        }
        numrecords++;

        ObjMgrFreeByEntityID (entityID);
      }

    } else {

      AsnReadVal (aip, atp, NULL);
    }
  }

  if (aop != NULL) {
    AsnCloseStruct (aop, cfp->atp_bsss, (Pointer) &(cfp->bss.seq_set));
    AsnCloseStruct (aop, cfp->bssp_atp, (Pointer) &(cfp->bss));
    AsnIoClose (aop);
  }

  AsnIoFree (aip, FALSE);

#ifdef OS_UNIX
  if (usedPopen) {
    pclose (fp);
  } else {
    FileClose (fp);
  }
#else
  FileClose (fp);
#endif
  if (cfp->logfp != NULL) {
    if (StringChr (cfp->report, 'c') == NULL) {
      fprintf (cfp->logfp, "\nTotal number of records %ld\n", (long) numrecords);
      if (StringDoesHaveText (longest)) {
        fprintf (cfp->logfp, "Longest processing time %ld seconds on %s\n",
                 (long) worsttime, longest);
      }
      fprintf (cfp->logfp, "Counts ");
      fprintf (cfp->logfp, "- %9ld RECS", (long) cfp->rawcounts.recs);
      fprintf (cfp->logfp, ", %9ld NUCS", (long) cfp->rawcounts.nucs);
      fprintf (cfp->logfp, ", %9ld PRTS", (long) cfp->rawcounts.prts);
      fprintf (cfp->logfp, ", %9ld OKAY", (long) cfp->rawcounts.okay);
      fprintf (cfp->logfp, ", %9ld NORM", (long) cfp->rawcounts.norm);
      fprintf (cfp->logfp, ", %9ld CLNR", (long) cfp->rawcounts.clnr);
      fprintf (cfp->logfp, ", %9ld OTHR", (long) cfp->rawcounts.othr);
      fprintf (cfp->logfp, ", %9ld MODR", (long) cfp->rawcounts.modr);
      fprintf (cfp->logfp, ", %9ld SLOC", (long) cfp->rawcounts.sloc);
      fprintf (cfp->logfp, ", %9ld PUBL", (long) cfp->rawcounts.publ);
      fprintf (cfp->logfp, ", %9ld AUTH", (long) cfp->rawcounts.auth);
      fprintf (cfp->logfp, ", %9ld SORT", (long) cfp->rawcounts.sort);
      fprintf (cfp->logfp, ", %9ld BSEC", (long) cfp->rawcounts.bsec);
      fprintf (cfp->logfp, ", %9ld GBBK", (long) cfp->rawcounts.gbbk);
      fprintf (cfp->logfp, ", %9ld TITL", (long) cfp->rawcounts.titl);
      fprintf (cfp->logfp, ", %9ld PACK", (long) cfp->rawcounts.pack);
      fprintf (cfp->logfp, ", %9ld MOVE", (long) cfp->rawcounts.move);
      fprintf (cfp->logfp, ", %9ld SSEC", (long) cfp->rawcounts.ssec);
      fprintf (cfp->logfp, "\n");
      fflush (cfp->logfp);
    }
  }
}

static void CleanupOneRecord (
  CharPtr filename,
  Pointer userdata
)

{
  CleanFlagPtr  cfp;
  CharPtr       ptr;
  SumDataPtr    sdp;

  if (StringHasNoText (filename)) return;
  cfp = (CleanFlagPtr) userdata;
  if (cfp == NULL) return;

  MemSet ((Pointer) &(cfp->rawcounts), 0, sizeof (CountData));
  MemSet ((Pointer) &(cfp->dbsums), 0, sizeof (DbSumData));

  if (StringChr (cfp->sourcedb, 'y') != NULL) {
    ptr = StringRChr (filename, DIRDELIMCHR);
    if (ptr != NULL) {
      ptr++;
      if (StringStr (ptr, "gbcon") != NULL ||
          StringStr (ptr, "gbest") != NULL ||
          StringStr (ptr, "gbgss") != NULL ||
          StringStr (ptr, "gbhtg") != NULL ||
          StringStr (ptr, "gbpat") != NULL ||
          StringStr (ptr, "gbsts") != NULL) return;
    }
  }

  if (cfp->batch) {
    ptr = StringRChr (filename, DIRDELIMCHR);
    if (ptr != NULL) {
      ptr++;
      if (StringDoesHaveText (cfp->firstfile)) {
        if (StringICmp (cfp->firstfile, ptr) == 0) {
          cfp->foundfirst = TRUE;
        }
        if (! cfp->foundfirst) return;
      }

      if (StringDoesHaveText (cfp->lastfile)) {
        if (cfp->foundlast) return;
        if (StringICmp (cfp->lastfile, ptr) == 0) {
          cfp->foundlast = TRUE;
        }
      }
    }

    CleanupMultipleRecord (filename, cfp);
  } else {
    CleanupSingleRecord (filename, cfp);
  }

  if (cfp->logfp != NULL) {
    if (StringChr (cfp->report, 'c') != NULL) {
      ptr = StringRChr (filename, DIRDELIMCHR);
      if (ptr != NULL) {
        ptr++;
        fprintf (cfp->logfp, "%s", ptr);
      }
      sdp = &(cfp->dbsums.genbank);
      if (sdp != NULL) {
        fprintf (cfp->logfp, "\t%ld\t%ld\t%ld", (long) sdp->recs, (long) sdp->nucs, (long) sdp->prts);
      }
      sdp = &(cfp->dbsums.embl);
      if (sdp != NULL) {
        fprintf (cfp->logfp, "\t%ld\t%ld\t%ld", (long) sdp->recs, (long) sdp->nucs, (long) sdp->prts);
      }
      sdp = &(cfp->dbsums.ddbj);
      if (sdp != NULL) {
        fprintf (cfp->logfp, "\t%ld\t%ld\t%ld", (long) sdp->recs, (long) sdp->nucs, (long) sdp->prts);
      }
      sdp = &(cfp->dbsums.refseq);
      if (sdp != NULL) {
        fprintf (cfp->logfp, "\t%ld\t%ld\t%ld", (long) sdp->recs, (long) sdp->nucs, (long) sdp->prts);
      }
      sdp = &(cfp->dbsums.other);
      if (sdp != NULL) {
        fprintf (cfp->logfp, "\t%ld\t%ld\t%ld", (long) sdp->recs, (long) sdp->nucs, (long) sdp->prts);
      }
      fprintf (cfp->logfp, "\t%ld\t%ld\t%ld", (long) cfp->rawcounts.recs, (long) cfp->rawcounts.nucs, (long) cfp->rawcounts.prts);
      fprintf (cfp->logfp, "\n");
      fflush (cfp->logfp);
    }
  }
}

static Boolean IsAllDigits (CharPtr str)

{
  CharPtr cp;

  if (StringHasNoText (str)) return FALSE;

  cp = str;
  while (*cp != 0 && isdigit (*cp)) {
    cp++;
  }
  if (*cp == 0) {
    return TRUE;
  } else {
    return FALSE;
  }
}

static Boolean PvtStrToLong (CharPtr str, Int4Ptr longval)

{
  Char     ch;
  Int2     i;
  Int2     len;
  Char     local [64];
  Boolean  nodigits;
  Boolean  rsult;
  long int     val;

  rsult = FALSE;
  if (longval != NULL) {
    *longval = (Int4) 0;
  }
  len = (Int2) StringLen (str);
  if (len != 0) {
    rsult = TRUE;
    nodigits = TRUE;
    for (i = 0; i < len; i++) {
      ch = str [i];
      if (ch == ' ' || ch == '+' || ch == '-') {
      } else if (ch < '0' || ch > '9') {
        rsult = FALSE;
      } else {
        nodigits = FALSE;
      }
    }
    if (nodigits) {
      rsult = FALSE;
    }
    if (rsult && longval != NULL) {
      StringNCpy_0 (local, str, sizeof (local));
      if (sscanf (local, "%ld", &val) == 1) {
        *longval = val;
      }
    }
  }
  return rsult;
}

static Int4 AccessionToGi (CharPtr string)

{
   Int4      gi;
   SeqIdPtr  sip;

   sip = SeqIdFromAccessionDotVersion (string);
   if (sip == NULL) return 0;
   gi = GetGIForSeqId (sip);
   SeqIdFree (sip);
   return gi;
}

static void ProcessAccessionList (
  CharPtr accnfile,
  CleanFlagPtr cfp
)

{
  AsnIoPtr     aop = NULL;
  BioseqPtr    bsp;
  Uint2        entityID;
  FileCache    fc;
  FILE         *fp;
  Char         line [1023];
  Char         path [PATH_MAX];
  SeqEntryPtr  sep;
  ErrSev       sev;
  SeqIdPtr     sip;
  CharPtr      str;
  Char         tmp [64];

  if (StringHasNoText (accnfile) || cfp == NULL) return;

  fp = FileOpen (accnfile, "r");
  if (fp == NULL) return;

  if (FileCacheSetup (&fc, fp)) {
    str = FileCacheGetString (&fc, line, sizeof (line));
    while (str != NULL) {
      TrimSpacesAroundString (str);

      bsp = NULL;
      entityID = 0;
      sep = NULL;
      sip = NULL;

      if (IsAllDigits (str)) {
        sprintf (tmp, "gi|%s", str);
        sip = SeqIdParse (tmp);
      } else {
        sip = SeqIdFromAccessionDotVersion (str);
      }

      if (sip != NULL) {
        sev = ErrSetMessageLevel (SEV_MAX);
        bsp = BioseqLockById (sip);
        ErrSetMessageLevel (sev);
        sip = SeqIdFree (sip);
      }

      if (bsp != NULL) {
        entityID = ObjMgrGetEntityIDForPointer (bsp);
        sep = GetTopSeqEntryForEntityID (entityID);

        if (sep != NULL) {

          path [0] = '\0';
          if (StringDoesHaveText (cfp->outfile)) {

            StringNCpy_0 (path, cfp->outfile, sizeof (path));

          } else if (StringDoesHaveText (cfp->results)) {

            StringNCpy_0 (path, cfp->results, sizeof (path));
            FileBuildPath (path, NULL, str);
            StringCat (path, ".asn");
          }

          entityID = ObjMgrGetEntityIDForChoice (sep);
          if (entityID > 0) {

            if (StringHasNoText (cfp->report) && StringDoesHaveText (path)) {
              aop = AsnIoOpen (path, "w");
            }

            DoCleanup (sep, entityID, cfp, aop, NULL, NULL);

            if (aop != NULL) {
              AsnIoFlush (aop);
              AsnIoClose (aop);
            }
          }
        }

        BioseqUnlock (bsp);
        SeqEntryFree (sep);
      }

      str = FileCacheGetString (&fc, line, sizeof (line));
    }
  }

  FileClose (fp);
}

/* Args structure contains command-line arguments */

typedef enum {
  p_argInputPath = 0,
  r_argOutputPath,
  i_argInputFile,
  o_argOutputFile,
  f_argFilter,
  x_argSuffix,
  j_argFirstFile,
  k_argLastFile,
  d_argSourceDb,
  a_argType,
  b_argBinary,
  c_argCompressed,
  L_argLogFile,
  R_argRemote,
  Q_argReport,
  S_argSelective,
  m_argFfMode,
  q_argFfDiff,
  n_argAsn2Flat,
  v_argAsnVal,
  K_argClean,
  U_argModernize,
  N_argLink,
  F_argFeat,
  D_argDesc,
  X_argMods,
  M_argMacro,
  T_argTaxonLookup,
  P_argPubLookup,
  A_argAccnFile,
#ifdef INTERNAL_NCBI_CLEANASN
  H_argHup,
#endif
} Arguments;

Args myargs [] = {
  {"Path to Files", NULL, NULL, NULL,
    TRUE, 'p', ARG_STRING, 0.0, 0, NULL},
  {"Path for Results", NULL, NULL, NULL,
    TRUE, 'r', ARG_STRING, 0.0, 0, NULL},
  {"Single Input File", "stdin", NULL, NULL,
    TRUE, 'i', ARG_FILE_IN, 0.0, 0, NULL},
  {"Single Output File", "stdout", NULL, NULL,
    TRUE, 'o', ARG_FILE_OUT, 0.0, 0, NULL},
  {"Substring Filter", NULL, NULL, NULL,
    TRUE, 'f', ARG_STRING, 0.0, 0, NULL},
  {"File Selection Suffix", ".ent", NULL, NULL,
    TRUE, 'x', ARG_STRING, 0.0, 0, NULL},
  {"First File Name", NULL, NULL, NULL,
    TRUE, 'j', ARG_STRING, 0.0, 0, NULL},
  {"Last File Name", NULL, NULL, NULL,
    TRUE, 'k', ARG_STRING, 0.0, 0, NULL},
  {"Source Database\n"
   "      a Any\n"
   "      g GenBank\n"
   "      e EMBL\n"
   "      d DDBJ\n"
   "      r RefSeq\n"
   "      n NCBI\n"
   "      x Exclude EMBL/DDBJ\n"
   "      y Exclude gbcon, gbest, gbgss, gbhtg, gbpat, gbsts\n", "a", NULL, NULL,
    TRUE, 'd', ARG_STRING, 0.0, 0, NULL},
  {"ASN.1 Type\n"
   "      a Any\n"
   "      e Seq-entry\n"
   "      b Bioseq\n"
   "      s Bioseq-set\n"
   "      m Seq-submit\n"
   "      t Batch Processing\n", "a", NULL, NULL,
    TRUE, 'a', ARG_STRING, 0.0, 0, NULL},
  {"Bioseq-set is Binary", "F", NULL, NULL,
    TRUE, 'b', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Bioseq-set is Compressed", "F", NULL, NULL,
    TRUE, 'c', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Log File", NULL, NULL, NULL,
    TRUE, 'L', ARG_FILE_OUT, 0.0, 0, NULL},
  {"Remote Fetching from ID", "F", NULL, NULL,
    TRUE, 'R', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Report\n"
   "      c Record Count\n"
   "      r ASN.1 BSEC Report\n"
   "      s ASN.1 SSEC Report\n"
   "      n NORM vs. SSEC Report\n"
   "      e PopPhyMutEco AutoDef Report\n"
   "      o Overlap Report\n"
   "      d Log SSEC Differences\n"
   "      g GenBank SSEC Diff\n"
   "      f asn2gb/asn2flat Diff\n"
   "      v Validator SSEC Diff\n"
   "      m Modernize Gene/RNA/PCR\n"
   "      u Unpublished Pub Lookup\n"
   "      p Published Pub Lookup\n"
   "      x Custom Scan\n", NULL, NULL, NULL,
    TRUE, 'Q', ARG_STRING, 0.0, 0, NULL},
  {"Selective Difference Filter\n"
   "      s SSEC\n"
   "      b BSEC\n"
   "      a Author\n"
   "      p Publication\n"
   "      l Location\n"
   "      r RNA\n"
   "      q Qualifier Sort Order\n"
   "      g Genbank Block\n"
   "      k Package CdRegion or Parts Features\n"
   "      m Move Publication\n"
   "      o Leave Duplicate Bioseq Publication\n"
   "      d Automatic Definition Line\n"
   "      e Pop/Phy/Mut/Eco Set Definition Line\n"
   "      (Capital Letters Skip)\n", NULL, NULL, NULL,
    TRUE, 'S', ARG_STRING, 0.0, 0, NULL},
  {"Flatfile Mode\n"
   "      r Release\n"
   "      e Entrez\n"
   "      s Sequin\n"
   "      d Dump\n", NULL, NULL, NULL,
    TRUE, 'm', ARG_STRING, 0.0, 0, NULL},
  {"ffdiff Executable", "/netopt/genbank/subtool/bin/ffdiff", NULL, NULL,
    TRUE, 'q', ARG_FILE_IN, 0.0, 0, NULL},
  {"asn2flat Executable", "/netopt/ncbi_tools/bin/asn2flat", NULL, NULL,
    TRUE, 'n', ARG_FILE_IN, 0.0, 0, NULL},
  {"asnval Executable", "/netopt/ncbi_tools/bin/asnval", NULL, NULL,
    TRUE, 'v', ARG_FILE_IN, 0.0, 0, NULL},
  {"Cleanup\n"
   "      b BasicSeqEntryCleanup\n"
   "      p C++ BasicCleanup\n"
   "      s SeriousSeqEntryCleanup\n"
   "      g GpipeSeqEntryCleanup\n"
   "      n Normalize Descriptor Order\n"
   "      u Remove NcbiCleanup User Objects\n", NULL, NULL, NULL,
    TRUE, 'K', ARG_STRING, 0.0, 0, NULL},
  {"Modernize\n"
   "      g Gene\n"
   "      r RNA\n"
   "      p PCR Primers\n", NULL, NULL, NULL,
    TRUE, 'U', ARG_STRING, 0.0, 0, NULL},
  {"Link\n"
   "      o LinkCDSmRNAbyOverlap\n"
   "      p LinkCDSmRNAbyProduct\n"
   "      r ReassignFeatureIDs\n"
   "      c ClearFeatureIDs\n", NULL, NULL, NULL,
    TRUE, 'N', ARG_STRING, 0.0, 0, NULL},
  {"Feature\n"
   "      u Remove User Object\n"
   "      d Remove db_xref\n"
   "      r Remove Redundant Gene xref\n", NULL, NULL, NULL,
    TRUE, 'F', ARG_STRING, 0.0, 0, NULL},
  {"Descriptor\n"
   "      t Remove Title\n"
   "      n Remove Nuc-Prot Set Title\n"
   "      e Remove Pop/Phy/Mut/Eco Set Title\n", NULL, NULL, NULL,
    TRUE, 'D', ARG_STRING, 0.0, 0, NULL},
  {"Miscellaneous\n"
   "      d Automatic Definition Line\n"
   "      e Pop/Phy/Mut/Eco Set Definition Line\n", NULL, NULL, NULL,
    TRUE, 'X', ARG_STRING, 0.0, 0, NULL},
  {"Macro File", NULL, NULL, NULL,
    TRUE, 'M', ARG_FILE_IN, 0.0, 0, NULL},
  {"Taxonomy Lookup", "F", NULL, NULL,
    TRUE, 'T', ARG_BOOLEAN, 0.0, 0, NULL},
  {"Publication Lookup\n"
   "      u Update PMID-only Publication\n"
   "      # Replace Unpublished With PMID\n", NULL, NULL, NULL,
    TRUE, 'P', ARG_STRING, 0.0, 0, NULL},
  {"Accession List File", NULL, NULL, NULL,
    TRUE, 'A', ARG_FILE_IN, 0.0, 0, NULL},
#ifdef INTERNAL_NCBI_CLEANASN
  {"Internal Access to HUP", "F", NULL, NULL,
    TRUE, 'H', ARG_BOOLEAN, 0.0, 0, NULL},
#endif
};

Int2 Main (void)

{
  ValNodePtr     action_list, parflat_list, vnp;
  AsnIoPtr       aip;
  Char           app [64], mode, type;
  CleanFlagData  cfd;
  CharPtr        directory, filter, accnfile, infile, logfile, outfile,
                 macro_file, results, str, suffix;
  Boolean        remote;
  time_t         runtime, starttime, stoptime;
#ifdef INTERNAL_NCBI_CLEANASN
  Boolean        hup;
#endif

  /* standard setup */

  ErrSetFatalLevel (SEV_MAX);
  ErrClearOptFlags (EO_SHOW_USERSTR);
  UseLocalAsnloadDataAndErrMsg ();
  ErrPathReset ();

  /* finish resolving internal connections in ASN.1 parse tables */

  if (! AllObjLoad ()) {
    Message (MSG_FATAL, "AllObjLoad failed");
    return 1;
  }
  if (! SubmitAsnLoad ()) {
    Message (MSG_FATAL, "SubmitAsnLoad failed");
    return 1;
  }
  if (! FeatDefSetLoad ()) {
    Message (MSG_FATAL, "FeatDefSetLoad failed");
    return 1;
  }
  if (! SeqCodeSetLoad ()) {
    Message (MSG_FATAL, "SeqCodeSetLoad failed");
    return 1;
  }
  if (! GeneticCodeTableLoad ()) {
    Message (MSG_FATAL, "GeneticCodeTableLoad failed");
    return 1;
  }

  parflat_list = Validate_ParFlat_GBFeat ();
  if (parflat_list != NULL) {
    Message (MSG_POSTERR, "Validate_ParFlat_GBFeat warnings");
    for (vnp = parflat_list; vnp != NULL; vnp = vnp->next) {
      str = (CharPtr) vnp->data.ptrvalue;
      if (StringHasNoText (str)) continue;
      Message (MSG_POSTERR, "%s", str);
    }
    ValNodeFreeData (parflat_list);
  }

  /* process command line arguments */

  sprintf (app, "cleanasn %s", CLEANASN_APPLICATION);
  if (! GetArgs (app, sizeof (myargs) / sizeof (Args), myargs)) {
    return 0;
  }

  MemSet ((Pointer) &cfd, 0, sizeof (CleanFlagData));

  directory = (CharPtr) myargs [p_argInputPath].strvalue;
  results = (CharPtr) myargs [r_argOutputPath].strvalue;
  if (StringHasNoText (results)) {
    results = directory;
  }
  accnfile = (CharPtr) myargs [A_argAccnFile].strvalue;
  infile = (CharPtr) myargs [i_argInputFile].strvalue;
  outfile = (CharPtr) myargs [o_argOutputFile].strvalue;
  filter = (CharPtr) myargs [f_argFilter].strvalue;
  suffix = (CharPtr) myargs [x_argSuffix].strvalue;

  cfd.batch = FALSE;
  cfd.binary = (Boolean) myargs [b_argBinary].intvalue;
  cfd.compressed = (Boolean) myargs [c_argCompressed].intvalue;
  cfd.type = 1;

  cfd.foundfirst = FALSE;
  cfd.foundlast = FALSE;
  cfd.sourcedb = myargs [d_argSourceDb].strvalue;

  str = myargs [a_argType].strvalue;
  TrimSpacesAroundString (str);
  if (StringDoesHaveText (str)) {
    type = str [0];
  } else {
    type = 'a';
  }

  type = TO_LOWER (type);
  switch (type) {
    case 'a' :
      cfd.type = 1;
      break;
    case 'e' :
      cfd.type = 2;
      break;
    case 'b' :
      cfd.type = 3;
      break;
    case 's' :
      cfd.type = 4;
      break;
    case 'm' :
      cfd.type = 5;
      break;
    case 't' :
      cfd.type = 1;
      cfd.batch = TRUE;
      break;
    default :
      cfd.type = 1;
      break;
  }

  remote = (Boolean) myargs [R_argRemote].intvalue;
  if (StringDoesHaveText (accnfile)) {
    remote = TRUE;
  }
#ifdef INTERNAL_NCBI_CLEANASN
  hup = (Boolean) myargs [H_argHup].intvalue;
#endif

  cfd.report = myargs [Q_argReport].strvalue;
  cfd.selective = myargs [S_argSelective].strvalue;
  cfd.ffdiff = myargs [q_argFfDiff].strvalue;
  cfd.asn2flat = myargs [n_argAsn2Flat].strvalue;
  cfd.asnval = myargs [v_argAsnVal].strvalue;

  str = myargs [m_argFfMode].strvalue;
  TrimSpacesAroundString (str);
  if (StringDoesHaveText (str)) {
    mode = str [0];
  } else {
    mode = 'e';
  }

  mode = TO_LOWER (mode);
  switch (mode) {
    case 'r' :
      cfd.ffmode = RELEASE_MODE;
      break;
    case 'e' :
      cfd.ffmode = ENTREZ_MODE;
      break;
    case 's' :
      cfd.ffmode = SEQUIN_MODE;
      break;
    case 'd' :
      cfd.ffmode = DUMP_MODE;
      break;
    default :
      cfd.ffmode = ENTREZ_MODE;
      break;
  }

  cfd.clean = myargs [K_argClean].strvalue;
  cfd.modernize = myargs [U_argModernize].strvalue;
  cfd.link = myargs [N_argLink].strvalue;
  cfd.feat = myargs [F_argFeat].strvalue;
  cfd.desc = myargs [D_argDesc].strvalue;
  cfd.mods = myargs [X_argMods].strvalue;
  cfd.taxon = (Boolean) myargs [T_argTaxonLookup].intvalue;
  cfd.pub = myargs [P_argPubLookup].strvalue;

  macro_file = myargs [M_argMacro].strvalue;
  if (StringDoesHaveText (macro_file)) {
    aip = AsnIoOpen (macro_file, "r");
    if (aip == NULL) {
      Message (MSG_FATAL, "Unable to open macro file '%s'", macro_file);
      return 1;
    }
    action_list = MacroActionListAsnRead (aip, NULL);
    AsnIoClose (aip);
    if (action_list == NULL) {
      Message (MSG_FATAL, "Unable to read macro file '%s'", macro_file);
    }
    cfd.action_list = action_list;
  }

  cfd.amp = AsnAllModPtr ();
  cfd.atp_bss = AsnFind ("Bioseq-set");
  cfd.atp_bsss = AsnFind ("Bioseq-set.seq-set");
  cfd.atp_se = AsnFind ("Bioseq-set.seq-set.E");
  cfd.atp_bsc = AsnFind ("Bioseq-set.class");
  cfd.bssp_atp = AsnLinkType (NULL, cfd.atp_bss);

  logfile = (CharPtr) myargs [L_argLogFile].strvalue;
  if (StringDoesHaveText (logfile)) {
    cfd.logfp = FileOpen (logfile, "w");
  }

  if (remote) {
#ifdef INTERNAL_NCBI_CLEANASN
    if (! PUBSEQBioseqFetchEnable ("cleanasn", FALSE)) {
      Message (MSG_POSTERR, "PUBSEQBioseqFetchEnable failed");
      return 1;
    }
    if (hup) {
      SmartFetchEnable ();
      TPASmartFetchEnable ();
      HUPFetchEnable ();
    }
#else
    PubSeqFetchEnable ();
#endif
  }

  if (remote || cfd.pub) {
    PubMedFetchEnable ();
  }

  /*
  if (cfd.logfp != NULL && StringChr (cfd.report, 'c') != NULL) {
    fprintf (cfd.logfp, "FILE\t\tGENBANK\t\t\tEMBL\t\t\tDDBJ\t\t\tREFSEQ\t\t\tOTHER\n");
    fprintf (cfd.logfp, "\tREC\tNUC\tPRT\tREC\tNUC\tPRT\tREC\tNUC\tPRT\tREC\tNUC\tPRT\tREC\tNUC\tPRT\n");
    fflush (cfd.logfp);
  }
  */

  starttime = GetSecs ();

  if (StringDoesHaveText (directory)) {
    if (StringHasNoText (cfd.report) && StringCmp (directory, results) == 0) {
      Message (MSG_POSTERR, "-r results path must be different than -p data path");
      if (cfd.logfp != NULL) {
        fprintf (cfd.logfp, "-r results path must be different than -p data path\n");
      }
    } else {

      cfd.firstfile = (CharPtr) myargs [j_argFirstFile].strvalue;
      cfd.lastfile = (CharPtr) myargs [k_argLastFile].strvalue;

      cfd.results = results;

      DirExplore (directory, filter, suffix, FALSE, CleanupOneRecord, (Pointer) &cfd);
    }

  } else if (StringDoesHaveText (accnfile)) {
    if (StringHasNoText (cfd.report) && StringHasNoText (results)) {
      Message (MSG_POSTERR, "-r results path must be set when -A accession list is used");
      if (cfd.logfp != NULL) {
        fprintf (cfd.logfp, "-r results path must be set when -A accession list is used\n");
      }
    } else {

      cfd.results = results;

      ProcessAccessionList (accnfile, (Pointer) &cfd);
    }

  } else if (StringDoesHaveText (infile) && StringDoesHaveText (outfile)) {

    cfd.outfile = outfile;

    CleanupOneRecord (infile, (Pointer) &cfd);
  }

  stoptime = GetSecs ();
  runtime = stoptime - starttime;

  if (cfd.logfp != NULL) {
    if (StringChr (cfd.report, 'c') == NULL) {
      fprintf (cfd.logfp, "\nFinished in %ld seconds\n", (long) runtime);
      fprintf (cfd.logfp, "Cumulative counts ");
      fprintf (cfd.logfp, "- %9ld RECS", (long) cfd.cumcounts.recs);
      fprintf (cfd.logfp, ", %9ld NUCS", (long) cfd.cumcounts.nucs);
      fprintf (cfd.logfp, ", %9ld PRTS", (long) cfd.cumcounts.prts);
      fprintf (cfd.logfp, ", %9ld OKAY", (long) cfd.cumcounts.okay);
      fprintf (cfd.logfp, ", %9ld NORM", (long) cfd.cumcounts.norm);
      fprintf (cfd.logfp, ", %9ld CLNR", (long) cfd.cumcounts.clnr);
      fprintf (cfd.logfp, ", %9ld OTHR", (long) cfd.cumcounts.othr);
      fprintf (cfd.logfp, ", %9ld MODR", (long) cfd.cumcounts.modr);
      fprintf (cfd.logfp, ", %9ld SLOC", (long) cfd.cumcounts.sloc);
      fprintf (cfd.logfp, ", %9ld PUBL", (long) cfd.cumcounts.publ);
      fprintf (cfd.logfp, ", %9ld AUTH", (long) cfd.cumcounts.auth);
      fprintf (cfd.logfp, ", %9ld SORT", (long) cfd.cumcounts.sort);
      fprintf (cfd.logfp, ", %9ld BSEC", (long) cfd.cumcounts.bsec);
      fprintf (cfd.logfp, ", %9ld GBBK", (long) cfd.cumcounts.gbbk);
      fprintf (cfd.logfp, ", %9ld TITL", (long) cfd.cumcounts.titl);
      fprintf (cfd.logfp, ", %9ld PACK", (long) cfd.cumcounts.pack);
      fprintf (cfd.logfp, ", %9ld MOVE", (long) cfd.cumcounts.move);
      fprintf (cfd.logfp, ", %9ld SSEC", (long) cfd.cumcounts.ssec);
      fprintf (cfd.logfp, "\n");
      fflush (cfd.logfp);
    }
    FileClose (cfd.logfp);
  }

  if (remote || cfd.pub) {
    PubMedFetchDisable ();
  }

  if (remote) {
#ifdef INTERNAL_NCBI_CLEANASN
    PUBSEQBioseqFetchDisable ();
#else
    PubSeqFetchDisable ();
#endif
  }

  return 0;
}

