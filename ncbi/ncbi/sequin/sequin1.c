/*   sequin1.c
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
* File Name:  sequin1.c
*
* Author:  Jonathan Kans
*
* Version Creation Date:   1/22/95
*
* $Revision: 6.729 $
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

#ifndef CODECENTER
static char *date_of_compilation = __DATE__;
static char *time_of_compilation = __TIME__;
#else
static char *date_of_compilation = "today";
static char *time_of_compilation = "now";
#endif

#include "sequin.h"
#include <vsm.h>
#include <vsmutil.h>
#include <valid.h>
#include <fstyle.h>
#include <biosrc.h>
#include <seqsub.h>
#include <cdrgn.h>
#include <import.h>
#include <medview.h>
#include <bspview.h>
#include <pubdesc.h>
#include <toasn3.h>
#include <utilpub.h>
#include <tofasta.h>
#include <saledit.h>
#include <salstruc.h>
#include <salfiles.h>
#include <salign.h>
#include <salsap.h>
#include <salutil.h>
#include <salpedit.h>
#include <salpanel.h>
#include <salptool.h>
#include <pobutil.h>
#include <accutils.h>
#include <netcnfg.h>
#include <objproj.h>
#include <suggslp.h>
#include <subutil.h>
#include <explore.h>
#include <actutils.h>
#include <pmfapi.h>
#include <accid1.h>
#include <ddvopen.h>
#include <dotseq.h>
#include <ingenwin.h>
#include <util/creaders/alnread.h>
#include <sqnutils.h>
#include <tax3api.h>
#include <validerr.h>
#include <algo/blast/api/blast_api.h>
#include <findrepl.h>

/* USE_SMARTNET */
#ifdef USE_SMARTNET
#include <smartnet.h>
#endif

/* USE_ENTREZ */
#include <accentr.h>

#include <objmime.h>
#include <mmdbapi.h>
/*
#ifndef WIN16
#include <cn3dentr.h>
#include <cn3dopen.h>
#endif
*/

#include <entrez.h>

/* USE_LOCAL */
#include <lsqfetch.h>

/* USE_MEDARCH */
#include <medarch.h>
#include <medutil.h>

#include <mla2api.h>


#ifdef USE_SPELL
#include <spellapi.h>
#endif

#ifdef OS_MAC
#include <Gestalt.h>
#endif

#define SEQ_APP_VER "10.05"

CharPtr SEQUIN_APPLICATION = SEQ_APP_VER;
CharPtr SEQUIN_SERVICES = NULL;
CharPtr SEQUIN_VERSION = NULL;

Boolean  useDesktop = FALSE;
Boolean  useEntrez = FALSE;
Boolean  useSeqFetch = FALSE;
Boolean  useIdLookup = FALSE;
Boolean  useLocal = FALSE;
Boolean  useBlast = FALSE;
Boolean  useMedarch = FALSE;
Boolean  newMedarch = FALSE;
Boolean  useTaxon = FALSE;
Boolean  allowDownload = FALSE;
Boolean  extraServices = FALSE;
Boolean  indexerVersion = FALSE;

CharPtr  genomeCenter = NULL;

Boolean  leaveAsOldAsn = FALSE;
Boolean  newAlignReader = TRUE;

SeqEntryPtr     globalsep = NULL;
Uint2           globalEntityID = 0;
Char            globalPath [PATH_MAX];

static SequinBlockPtr  globalsbp = NULL;

static Boolean useOldGraphicView = FALSE;
static Boolean useOldAlignmentView = FALSE;
static Boolean useOldSequenceView = FALSE;
/* static Boolean useUdv = FALSE; */

static Boolean gphviewscorealigns = FALSE;

ForM  startupForm = NULL;

#ifdef WIN_MAC
Boolean  termListUp = FALSE;
Boolean  docSumUp = FALSE;
Boolean  bioseqViewUp = FALSE;
#endif

#ifdef WIN_MAC
static IteM  openItem = NULL;
static IteM  closeItem = NULL;
static IteM  importItem = NULL;
static IteM  exportItem = NULL;
static IteM  duplicateViewItem = NULL;
static IteM  saveItem = NULL;
static IteM  saveAsItem = NULL;
static IteM  restoreItem = NULL;
static IteM  prepareItem = NULL;
static IteM  submitItem = NULL;
static IteM  loadUidItem = NULL;
static IteM  saveUidItem = NULL;
static IteM  printItem = NULL;

static IteM  undoItem = NULL;
static IteM  cutItem = NULL;
static IteM  copyItem = NULL;
static IteM  pasteItem = NULL;
static IteM  deleteItem = NULL;
static IteM  duplicateItem = NULL;

static IteM  orfItem = NULL;
static IteM  aluItem = NULL;
static IteM  targetItem = NULL;
static IteM  findItem = NULL;
static IteM  findFFItem = NULL;
static IteM  findGeneItem = NULL;
static IteM  findProtItem = NULL;
static IteM  findPosItem = NULL;
static IteM  validateItem = NULL;
static MenU  validateMenu = NULL;
static MenU  vecscreenMenu = NULL;
static IteM  spellItem = NULL;
static IteM  vectorScreenItem = NULL;
static IteM  cddBlastItem = NULL;
static MenU  cddSearchMenu = NULL;
static IteM  cddSearchItem = NULL;
static IteM  editsequenceitem = NULL;
static IteM  editseqalignitem = NULL;
static IteM  editseqsubitem = NULL;
static IteM  edithistoryitem = NULL;
static MenU  updateSeqMenu = NULL;
static MenU  updateSeqMenuIndexer = NULL;
static MenU  extendSeqMenu = NULL;
static MenU  addSeqMenu = NULL;
static IteM  featPropItem = NULL;
static IteM  parseFileItem = NULL;
static IteM  updalignitem = NULL;
static IteM  updalignidxitem = NULL;

static IteM  docsumfontItem = NULL;
static IteM  displayfontItem = NULL;
static IteM  preferencesItem = NULL;
static IteM  clearUnusedItem = NULL;
static IteM  legendItem = NULL;
static ChoicE  queryChoice = NULL;
static ChoicE  neighborChoice = NULL;
static IteM  oldAsnItem = NULL;

extern IteM  addSecondaryItem;
IteM  addSecondaryItem = NULL;
#endif

static Int2  startupStyle = 0;

static ForM  termListForm = NULL;
static ForM  docSumForm = NULL;

static Boolean  loadSaveUidListOK = FALSE;

static MedlineViewProcs    medviewprocs;
SeqViewProcs        seqviewprocs;
static EntrezGlobals       entrezglobals;

static SeqEditViewProcs    seqedprocs;
static StdEditorProcs      stdedprocs;
static StdEditorProcs      valdtrprocs;
static TextViewProcs       txtviewprocs;
static PubdescEditProcs    pubedprocs;
static BioSourceEditProcs  biosrcedprocs;

/*
static PRGD  prgdDict = NULL;
*/

static Boolean  workbenchMode = FALSE;
static Boolean  subtoolMode = FALSE;
static Boolean  stdinMode = FALSE;
static Boolean  bioseqsetMode = FALSE;
static Boolean  binseqentryMode = FALSE;
static Boolean  entrezMode = FALSE;
static Boolean  nohelpMode = FALSE;
static Boolean  backupMode = FALSE;
static Uint2    subtoolDatatype = 0;
static Uint2    subtoolEntityID = 0;

static Boolean  smartnetMode = FALSE;
static Boolean  dirsubMode = FALSE;

#ifdef WIN_MAC
static MenU     newDescMenu = NULL;
static MenU     newFeatMenu = NULL;
static MenU     advTableMenu = NULL;
static IteM     sucItem = NULL;
static MenU     newPubMenu = NULL;
static MenU     batchApplyMenu = NULL;
static MenU     batchEditMenu = NULL;
static MenU     specialMenu = NULL;
static MenU     projectsMenu = NULL;
static MenU     analysisMenu = NULL;
static Boolean  initialFormsActive = FALSE;
#endif

static ForM  initSubmitForm = NULL;
static ForM  formatForm = NULL;

static Int2     subtoolTimerLimit = 100;
static Int2     subtoolTimerCount = 0;
static Boolean  subtoolRecordDirty = FALSE;

static Boolean  testLatLonSubregion = FALSE;
static Boolean  strictLatLonCountry = FALSE;


#ifdef USE_SMARTNET
static Int4 SMWriteBioseqObj(VoidPtr bio_data, SMUserDataPtr sm_usr_data, 
                             VoidPtr data);
static Int4 SMReadBioseqObj(VoidPtr data, CharPtr buffer, 
                            Int4 length, void* fd);
#define SMART_KEY 1313
#define DUMB_KEY 1314
#endif


static FormatBlock globalFormatBlock = {SEQ_PKG_SINGLE, SEQ_FMT_FASTA, 0, SEQ_ORIG_SUBMISSION};

ForM  helpForm = NULL;

static CharPtr validFailMsg =
"Submission failed validation test.  Continue?\n\
(Choose Validate in the Search menu to see errors.)";


extern Int2 GetSequinAppParam (CharPtr section, CharPtr type, CharPtr dflt, CharPtr buf, Int2 buflen)

{
  Int2  rsult;

  rsult = GetAppParam ("SEQUINCUSTOM", section, type, NULL, buf, buflen);
  if (rsult) return rsult;
  rsult = GetAppParam ("SEQUIN", section, type, dflt, buf, buflen);
  return rsult;
}

extern Boolean WriteSequinAppParam (CharPtr section, CharPtr type, CharPtr value)

{
  MsgAnswer  ans;

  if (indexerVersion) {
    if (section == NULL) {
      section = "";
    }
    if (type == NULL) {
      type = "";
    }
    if (value == NULL) {
      value = "";
    }
    ans = Message (MSG_OKC, "Writing to local .sequinrc configuration file - [%s] %s = %s", section, type, value);
    if (ans == ANS_CANCEL) return FALSE;
  }
  return SetAppParam ("SEQUIN", section, type, value);
}

static Boolean SetSequinAppParam (CharPtr section, CharPtr type, CharPtr value)

{
  Char  tmp [32];

  if (GetAppParam ("SEQUINCUSTOM", section, type, NULL, tmp, sizeof (tmp) - 1)) {
    return SetAppParam ("SEQUINCUSTOM", section, type, value);
  }
  return WriteSequinAppParam (section, type, value);
}

static void SetSequinAppParamTF (CharPtr section, CharPtr type, Boolean value)

{
  if (value) {
    SetSequinAppParam (section, type, "TRUE");
  } else {
    SetSequinAppParam (section, type, "FALSE");
  }
}

static void CheckForCookedBioseqs (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BoolPtr    bp;
  BioseqPtr  bsp;

  if (sep == NULL) return;
  if (! IS_Bioseq (sep)) return;
  bp = (BoolPtr) mydata;
  if (bp == NULL) return;
  bsp = (BioseqPtr) sep->data.ptrvalue;
  if (bsp == NULL) return;
  if (bsp->repr != Seq_repr_raw && bsp->repr != Seq_repr_seg &&
  	bsp->repr != Seq_repr_delta && bsp->repr != Seq_repr_virtual) {
    *bp = FALSE;
  }
}

static void TaxonValidate (SeqEntryPtr sep, ValidStructPtr vsp);


typedef enum {
  eOkToWriteEntity_Cancel = 0,
  eOkToWriteEntity_Continue,
  eOkToWriteEntity_Validate
} EOkToWriteEntity;

static EOkToWriteEntity GetValidationCancelContinue (ValidStructPtr vsp, Boolean allow_review)
{
  WindoW w;
  GrouP  h, c, prompts;
  PrompT p1, p2;
  ButtoN b;
  ModalAcceptCancelData acd;
  CharPtr msg, msg_format = "Reject %d, Error %d, Warning %d, Info %d";
  EOkToWriteEntity rval = eOkToWriteEntity_Cancel;
  
  acd.accepted = FALSE;
  acd.cancelled = FALSE;
  acd.third_option = FALSE;
  
  w = ModalWindow(-20, -13, -10, -10, NULL);
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  
  prompts = HiddenGroup (h, -1, 0, NULL);
  SetGroupSpacing (prompts, 10, 10);
  p1 = StaticPrompt (prompts, "Submission failed validation test with:", 0, 0, programFont, 'l');

  msg = (CharPtr) MemNew (sizeof (Char) * (StringLen (msg_format) + 75));
  sprintf (msg, msg_format, 
              (int) vsp->errors [4], (int) vsp->errors [3],
              (int) vsp->errors [2], (int) vsp->errors [1]);
  p2 = StaticPrompt (h, msg, 0, 0, programFont, 'l'); 
  msg = MemFree (msg);
  AlignObjects (ALIGN_CENTER, (HANDLE) p1, (HANDLE) p2, NULL);

  c = HiddenGroup (h, 3, 0, NULL);
  SetGroupSpacing (c, 10, 10);
  b = PushButton (c, "Continue", ModalAcceptButton);
  SetObjectExtra (b, &acd, NULL);
  if (allow_review) {
    b = PushButton (c, "Review Errors", ModalThirdOptionButton);
    SetObjectExtra (b, &acd, NULL);
  }
  b = PushButton (c, "Cancel", ModalCancelButton);
  SetObjectExtra (b, &acd, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) prompts, (HANDLE) c, NULL);
  
  Show(w); 
  Select (w);
  while (!acd.accepted && ! acd.cancelled && ! acd.third_option)
  {
    ProcessExternalEvent ();
    Update ();
  }
  ProcessAnEvent ();
  Remove (w);
  if (acd.third_option) 
  {
    rval = eOkToWriteEntity_Validate;
  } 
  else if (acd.accepted)
  {
    rval = eOkToWriteEntity_Continue;
  }
  else
  {
    rval = eOkToWriteEntity_Cancel;
  }
  return rval;
}


static EOkToWriteEntity OkayToWriteTheEntity (Uint2 entityID, ForM f, Boolean allow_review)

{
  Boolean          allRawOrSeg = TRUE;
  EOkToWriteEntity rval = eOkToWriteEntity_Continue;
  Int2             errors;
  Int2             j;
  ErrSev           oldErrSev;
  SeqEntryPtr      sep;
  Char             str [32];
  ValidStructPtr   vsp;

  if (entityID < 1) return 0;
  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep == NULL) return 0;

  if (!FixSpecialCharacters (entityID)) return 0;

  if (GetSequinAppParam ("PREFERENCES", "ASKBEFOREVALIDATE", NULL, str, sizeof (str))) {
    if (StringICmp (str, "TRUE") == 0) {
      if (! (subtoolMode ||smartnetMode || backupMode) ) {
        if (Message (MSG_YN, "Do you wish to validate this entry?") == ANS_NO) return 1;
      }
    }
  }
  WatchCursor ();
  Update ();
  vsp = ValidStructNew ();
  if (vsp != NULL) {
    /*SetChecklistValue (checklistForm, 6);*/
    SeqEntryExplore (sep, (Pointer) (&allRawOrSeg), CheckForCookedBioseqs);
    if (allRawOrSeg) {
      vsp->useSeqMgrIndexes = TRUE;
    }
    if (indexerVersion) {
      vsp->alwaysRequireIsoJTA = TRUE;
      if (smartnetMode) {
        vsp->farFetchCDSproducts = TRUE;
        vsp->farFetchMRNAproducts = TRUE;
      }
    }
    oldErrSev = ErrSetMessageLevel (SEV_MAX);
    vsp->validateAlignments = TRUE;
    vsp->alignFindRemoteBsp = TRUE;
    vsp->doSeqHistAssembly = FALSE;
    if (smartnetMode) {
      vsp->doSeqHistAssembly = TRUE;
      vsp->farIDsInAlignments = TRUE;
      if (useEntrez) {
        /*
        LookupFarSeqIDs (sep, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE);
        vsp->inferenceAccnCheck = TRUE;
        */
      }
      vsp->testLatLonSubregion = testLatLonSubregion;
      vsp->strictLatLonCountry = strictLatLonCountry;
    }
    for (j = 0; j < 6; j++) {
      vsp->errors [j] = 0;
    }
    ValidateSeqEntry (sep, vsp);
    if (indexerVersion && useEntrez) {
      TaxonValidate (sep, vsp);
    }
    ErrSetMessageLevel (oldErrSev);
    ErrClear ();
    ErrShow ();
    errors = 0;
    if (subtoolMode || smartnetMode || backupMode) {
      for (j = 0; j < 6; j++) {
        errors += vsp->errors [j];
      }
    } else {
      for (j = 3; j < 6; j++) {
        errors += vsp->errors [j];
      }
    }
    UseWindow ((WindoW) f);
    if (errors > 0) {
      ArrowCursor ();
      Update ();
      if (subtoolMode || smartnetMode || backupMode) {
        rval = GetValidationCancelContinue (vsp, allow_review);
      } else {
        if (Message (MSG_OKC, validFailMsg) == ANS_CANCEL) {
          rval = eOkToWriteEntity_Cancel;
        }
      }
    }
    ValidStructFree (vsp);
  }
  ArrowCursor ();
  Update ();
  return rval;
}

static void ReplaceString (CharPtr PNTR target, CharPtr newstr)

{
  if (target == NULL) return;
  MemFree (*target);
  *target = StringSaveNoNull (newstr);
}

static void UncompressBsps (BioseqPtr bsp, Pointer userdata)

{
  if (ISA_na (bsp->mol)) {
    BioseqConvert (bsp, Seq_code_iupacna);
  } else if (ISA_aa (bsp->mol)) {
    BioseqConvert (bsp, Seq_code_ncbieaa);
  }
}

extern Boolean WriteTheEntityID (Uint2 entityID, CharPtr path, Boolean binary)

{
  AsnIoPtr       aip;
  BioseqPtr      bsp;
  BioseqSetPtr   bssp;
  ObjMgrDataPtr  omdp;
  Boolean        rsult;
  SeqEntryPtr    sep;
  SeqSubmitPtr   ssp;
  Char           str [16];
#ifdef WIN_MAC
  FILE           *f;
#endif

  rsult = FALSE;
  if (entityID < 1) return rsult;
  if (path == NULL || path [0] == '\0') return rsult;
  ssp = NULL;
  sep = NULL;
  bsp = NULL;
  bssp = NULL;
  omdp = ObjMgrGetData (entityID);
  if (omdp == NULL) return rsult;
  WatchCursor ();
  Update ();
#ifdef WIN_MAC
  f = FileOpen (path, "r");
  if (f != NULL) {
    FileClose (f);
  } else {
    FileCreate (path, "TEXT", "ttxt");
  }
#endif
  if (GetSequinAppParam ("PREFERENCES", "UNCOMPRESS", NULL, str, sizeof (str))) {
    if (StringICmp (str, "TRUE") == 0) {
      sep = GetTopSeqEntryForEntityID (entityID);
      VisitBioseqsInSep (sep, NULL, UncompressBsps);
      /*
      SeqEntryConvert (sep, Seq_code_iupacna);
      SeqEntryConvert (sep, Seq_code_ncbieaa);
      */
    }
  }
  if (binseqentryMode) {
    aip = AsnIoOpen (path, "w");
    sep = GetTopSeqEntryForEntityID (entityID);
    rsult = SeqEntryAsnWrite (sep, aip, NULL);
    AsnIoClose (aip);
    ArrowCursor ();
    Update ();
    return rsult;
  }
  sep = NULL;
  if (binary) {
    aip = AsnIoOpen (path, "wb");
  } else {
    aip = AsnIoOpen (path, "w");
  }
  if (aip != NULL) {
    switch (omdp->datatype) {
      case OBJ_SEQSUB :
        ssp = (SeqSubmitPtr) omdp->dataptr;
        if (ssp != NULL && ssp->datatype == 1) {
          rsult = SeqSubmitAsnWrite (ssp, aip, NULL);
        }
        break;
      case OBJ_BIOSEQ :
        sep = (SeqEntryPtr) omdp->choice;
        if (sep != NULL && sep->choice == 1 && sep->data.ptrvalue != NULL) {
          if (subtoolMode || bioseqsetMode) {
            bsp = (BioseqPtr) sep->data.ptrvalue;
            rsult = BioseqAsnWrite (bsp, aip, NULL);
          } else {
            rsult = SeqEntryAsnWrite (sep, aip, NULL);
          }
        }
        break;
      case OBJ_BIOSEQSET :
        sep =  (SeqEntryPtr) omdp->choice;
        if (sep != NULL && sep->choice == 2 && sep->data.ptrvalue != NULL) {
          if (subtoolMode || bioseqsetMode) {
            bssp = (BioseqSetPtr) sep->data.ptrvalue;
            rsult = BioseqSetAsnWrite (bssp, aip, NULL);
          } else {
            rsult = SeqEntryAsnWrite (sep, aip, NULL);
          }
        }
        break;
      case OBJ_SEQENTRY :
        sep =  (SeqEntryPtr) omdp->choice;
        if (sep != NULL) {
          rsult = SeqEntryAsnWrite (sep, aip, NULL);
        }
        break;
      default :
        break;
    }
    AsnIoClose (aip);
    if (! smartnetMode) {
      ObjMgrSetDirtyFlag (entityID, FALSE);
    }
  }
  ArrowCursor ();
  Update ();
  return rsult;
}

static ValNodePtr ExtractGivenSeqDescrUserObject (ValNodePtr PNTR headptr, CharPtr str, CharPtr cls)

{
  Boolean        extract_it;
  ValNodePtr     last = NULL, vnp;
  ObjectIdPtr    oip;
  UserObjectPtr  uop;

  if (headptr == NULL) return NULL;
  vnp = *headptr;

  while (vnp != NULL) {
    extract_it = FALSE;
    if (vnp->choice == Seq_descr_user) {
      uop = (UserObjectPtr) vnp->data.ptrvalue;
      if (uop != NULL) {
        if (StringDoesHaveText (cls)) {
          if (StringICmp (uop->_class, cls) == 0) {
            extract_it = TRUE;
          }
        }
        if (StringDoesHaveText (str)) {
          oip = uop->type;
          if (oip != NULL) {
            if (StringICmp (oip->str, str) == 0) {
              extract_it = TRUE;
            }
          }
        }
      }
    }
    if (extract_it) {
      if (last == NULL) {
        *headptr = vnp->next;
      } else {
        last->next = vnp->next;
      }
      vnp->next = NULL;
      return vnp;
    } else {
      last = vnp;
      vnp = vnp->next;
    }
  }

  return NULL;
}

typedef struct propgenbankdata {
  Boolean  ask;
  Boolean  asked;
  Boolean  bail;
  Boolean  changed;
} PropGenbankData, PNTR PropGenBankPtr;

static void DoPropagateFromGenBankBioseqSet (
  BioseqSetPtr seqset,
  Pointer userdata
)

{
  BioseqPtr       bsp;
  BioseqSetPtr    bssp;
  PropGenBankPtr  pgp;
  SeqEntryPtr     seqentry;
  ValNodePtr      smartuserobj;
  ValNodePtr      sourcedescr;
  UserObjectPtr   uop;

  if (seqset == NULL) return;
  if (seqset->_class != BioseqseqSet_class_genbank) return;
  pgp = (PropGenBankPtr) userdata;
  if (pgp == NULL) return;

  seqentry = seqset->seq_set;
  sourcedescr = seqset->descr;
  if (sourcedescr == NULL) return;

  /* if only descriptor is tracking user object, skip */
  if (sourcedescr->next == NULL && sourcedescr->choice == Seq_descr_user) {
    uop = (UserObjectPtr) sourcedescr->data.ptrvalue;
    if (uop != NULL && StringICmp (uop->_class, "SMART_V1.0") == 0) return;
  }

  /* optionally ask if propagation is desired */
  if (pgp->ask) {
    if (! pgp->asked) {
      if (Message (MSG_YN, "Propagate descriptors from top-level set?") == ANS_NO) {
        pgp->bail = TRUE;
      }
      pgp->asked = TRUE;
    }
  }
  if (pgp->bail) return;

  /* disconnect descriptors from parent bssp */
  seqset->descr = NULL;

  /* extract tracking user object */
  smartuserobj = ExtractGivenSeqDescrUserObject (&sourcedescr, NULL, "SMART_V1.0");

  while (seqentry != NULL) {
    if (seqentry->data.ptrvalue != NULL) {
      if (seqentry->choice == 1) {
        bsp = (BioseqPtr) seqentry->data.ptrvalue;
        ValNodeLink (&(bsp->descr),
                     AsnIoMemCopy ((Pointer) sourcedescr,
                                   (AsnReadFunc) SeqDescrAsnRead,
                                   (AsnWriteFunc) SeqDescrAsnWrite));
      } else if (seqentry->choice == 2) {
        bssp = (BioseqSetPtr) seqentry->data.ptrvalue;
        ValNodeLink (&(bssp->descr),
                     AsnIoMemCopy ((Pointer) sourcedescr,
                                   (AsnReadFunc) SeqDescrAsnRead,
                                   (AsnWriteFunc) SeqDescrAsnWrite));
      }
      pgp->changed = TRUE;
    }
    seqentry = seqentry->next;
  }

  /* free extracted original descriptors now that copies are propagated */
  SeqDescrFree (sourcedescr);

  /* restore tracking user object */
  if (smartuserobj != NULL) {
    ValNodeLink (&(seqset->descr), smartuserobj);
  }

  /* recurse */
  VisitSetsInSet (seqset, userdata, DoPropagateFromGenBankBioseqSet);
}

extern Boolean PropagateFromGenBankBioseqSet (SeqEntryPtr sep, Boolean ask)

{
  BioseqSetPtr     bssp;
  PropGenbankData  pdp;

  if (sep == NULL) return FALSE;
  if (! IS_Bioseq_set (sep)) return FALSE;

  bssp = (BioseqSetPtr) sep->data.ptrvalue;
  if (bssp == NULL) return FALSE;
  if (bssp->_class != BioseqseqSet_class_genbank) return FALSE;

  MemSet ((Pointer) &pdp, 0, sizeof (PropGenbankData));
  pdp.ask = ask;
  pdp.asked = FALSE;
  pdp.bail = FALSE;
  pdp.changed = FALSE;

  DoPropagateFromGenBankBioseqSet (bssp, (Pointer) &pdp);

  return pdp.changed;
}

static void ForcePropagate (IteM i)

{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = (BaseFormPtr) currentFormDataPtr;
#else
  bfp = (BaseFormPtr) GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  PropagateFromGenBankBioseqSet (sep, FALSE);
  NormalizeDescriptorOrder (sep);
  ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
}

#define SEQUIN_EDIT_TEMP_FILE "sequinEdit.temp"
#define SEQUIN_EDIT_BACK_FILE "sequinEdit.back"
#define SEQUIN_EDIT_PREV_FILE "sequinEdit.prev"
#define SEQUIN_EDIT_ARCH_FILE "sequinEdit.arch"

static void SubtoolModeTimerProc (void)

{
  ObjMgrDataPtr  omdp;

  subtoolTimerCount++;
  if (subtoolTimerCount > subtoolTimerLimit) {
    subtoolTimerCount = 0;
    if (subtoolRecordDirty) {
      omdp = ObjMgrGetData (subtoolEntityID);
      if (omdp != NULL) {
        if (WriteTheEntityID (subtoolEntityID, SEQUIN_EDIT_TEMP_FILE, FALSE)) {
          FileRemove (SEQUIN_EDIT_PREV_FILE);
          FileRename (SEQUIN_EDIT_BACK_FILE, SEQUIN_EDIT_PREV_FILE);
          FileRename (SEQUIN_EDIT_TEMP_FILE, SEQUIN_EDIT_BACK_FILE);
        } else {
          Message (MSG_POSTERR, "Unable to save automatic temporary file");
        }
      }
      subtoolRecordDirty = FALSE;
    }
  }
  SequinCheckSocketsProc ();
}

static Int2 LIBCALLBACK SubtoolModeMsgFunc (OMMsgStructPtr ommsp)

{
  switch (ommsp->message) {
    case OM_MSG_DEL :
    case OM_MSG_CREATE :
    case OM_MSG_UPDATE :
      subtoolRecordDirty = TRUE;
      break;
    default :
      break;
  }
  return OM_MSG_RET_OK;
}

static void BackupModeTimerProc (void)

{
  ObjMgrDataPtr  omdp;

  subtoolTimerCount++;
  if (subtoolTimerCount > subtoolTimerLimit) {
    subtoolTimerCount = 0;
    if (subtoolRecordDirty && subtoolEntityID > 0) {
      omdp = ObjMgrGetData (subtoolEntityID);
      if (omdp != NULL) {
        if (WriteTheEntityID (subtoolEntityID, SEQUIN_EDIT_TEMP_FILE, FALSE)) {
          FileRemove (SEQUIN_EDIT_PREV_FILE);
          FileRename (SEQUIN_EDIT_BACK_FILE, SEQUIN_EDIT_PREV_FILE);
          FileRename (SEQUIN_EDIT_TEMP_FILE, SEQUIN_EDIT_BACK_FILE);
        } else {
          Message (MSG_POSTERR, "Unable to save automatic temporary file");
        }
      }
      subtoolRecordDirty = FALSE;
    }
  }
  SequinCheckSocketsProc ();
}

static Int2 LIBCALLBACK BackupModeMsgFunc (OMMsgStructPtr ommsp)

{
  switch (ommsp->message) {
    case OM_MSG_DEL :
    case OM_MSG_CREATE :
    case OM_MSG_UPDATE :
      subtoolRecordDirty = TRUE;
      break;
    default :
      break;
  }
  return OM_MSG_RET_OK;
}

static void GetDefaultTitleFromForm (ForM f, CharPtr str, size_t maxsize, CharPtr filepath)

{
  Char     ch;
  Char     dfault [64];
  Int2     j;
  Int2     k;
  CharPtr  ptr;

  if (f != NULL && str != NULL && maxsize > 0) {
    dfault [0] = '\0';
    if (StringHasNoText (filepath)) {
      GetTitle (f, dfault, sizeof (dfault));
    } else {
      ptr = StringRChr (filepath, DIRDELIMCHR);
      if (ptr != NULL) {
        ptr++;
        StringNCpy_0 (dfault, ptr, sizeof (dfault));
      } else {
        StringNCpy_0 (dfault, filepath, sizeof (dfault));
      }
    }
    j = 0;
    k = 0;
    ch = dfault [j];
    while (j < sizeof (dfault) && ch != '\0') {
      if (ch <= ' ') {
        j++;
      } else {
        dfault [k] = dfault [j];
        k++;
        j++;
      }
      ch = dfault [j];
    }
    dfault [k] = '\0';
#ifdef WIN_MSWIN
    j = 0;
    ch = dfault [j];
    while (j < sizeof (dfault) && ch != '\0') {
      if (ch == '_' || IS_ALPHANUM (ch)) {
        j++;
        ch = dfault [j];
      } else {
        ch = '\0';
      }
    }
    dfault [j] = '\0';
#endif
    StringNCpy_0 (str, dfault, maxsize);
  }
}

static CharPtr gbsub = "gb-sub@ncbi.nlm.nih.gov";
static CharPtr emblsub = "datasubs@ebi.ac.uk";
static CharPtr ddbjsub = "ddbjsub@ddbj.nig.ac.jp";

static CharPtr gbupd = "gb-admin@ncbi.nlm.nih.gov";
static CharPtr emblupd = "update@ebi.ac.uk";
static CharPtr ddbjupd = "ddbjupdt@ddbj.nig.ac.jp";

static CharPtr ReturnSubmissionEmailAddress (Uint2 entityID)

{
  ObjMgrDataPtr   omdp;
  CharPtr         rsult;
  SubmitBlockPtr  sbp;
  SeqSubmitPtr    ssp;
  Char            str [32];
  Boolean         update = FALSE;

  rsult = gbsub;
  omdp = ObjMgrGetData (entityID);
  if (omdp != NULL && omdp->datatype == OBJ_SEQSUB) {
    ssp = (SeqSubmitPtr) omdp->dataptr;
    if (ssp != NULL && ssp->datatype == 1) {
      sbp = ssp->sub;
      if (sbp != NULL && sbp->subtype == 2) {
        update = TRUE;
      }
    }
  }
  if (GetAppParam ("SEQUIN", "PREFERENCES", "DATABASE", NULL, str, sizeof (str))) {
    if (StringICmp (str, "GenBank") == 0) {
      if (update) {
        rsult = gbupd;
      } else {
        rsult = gbsub;
      }
    } else if (StringICmp (str, "EMBL") == 0) {
      if (update) {
        rsult = emblupd;
      } else {
        rsult = emblsub;
      }
    } else if (StringICmp (str, "DDBJ") == 0) {
      if (update) {
        rsult = ddbjupd;
      } else {
        rsult = ddbjsub;
      }
    }
  }
  return rsult;
}

static void DoRemoveAlignmentFromRecord (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  SeqAnnotPtr   nextsap;
  Pointer PNTR  prevsap;
  SeqAnnotPtr   sap;

  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sap = bsp->annot;
    prevsap = (Pointer PNTR) &(bsp->annot);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sap = bssp->annot;
    prevsap = (Pointer PNTR) &(bssp->annot);
  } else return;
  while (sap != NULL) {
    nextsap = sap->next;
    if (sap->type == 2) {
      *(prevsap) = sap->next;
      sap->next = NULL;
      SeqAnnotFree (sap);
    } else {
      prevsap = (Pointer PNTR) &(sap->next);
    }
    sap = nextsap;
  }
}

extern void SubmitToNCBI (IteM i);

static void MissingAnnotCallback (BioseqPtr bsp, Pointer userdata)
{
  BoolPtr p_missing;
  SeqFeatPtr sfp;
  SeqMgrFeatContext fcontext;

  if (bsp == NULL || ISA_aa (bsp->mol) || userdata == NULL) {
    return;
  }
  p_missing = (BoolPtr) userdata;
  if (*p_missing) {
    return;
  }
  sfp = SeqMgrGetNextFeature (bsp, NULL, 0, 0, &fcontext);
  if (sfp == NULL) {
    *p_missing = TRUE;
  }
}

  
static Boolean IsAnySequenceMissingAnnotation (SeqEntryPtr sep)
{
  Boolean rval = FALSE;

  VisitBioseqsInSep (sep, &rval, MissingAnnotCallback);
  return rval;
}


static void PrepareSeqSubmitProc (IteM i)

{
  BaseFormPtr  bfp;
  Char         dfault [64];
  Char         path [PATH_MAX];
  CharPtr      ptr;
  SeqEntryPtr  sep;
  CharPtr      str, email_address;
  Boolean      update;
  CharPtr      fmt_file = "Submission is now written.  Please e-mail '%s' to %s.%s";
  CharPtr      fmt_no_file = "Submission is now written.  Please e-mail to %s.%s";
  CharPtr      missing_annot = "  Please include a brief summary of your submission within your correspondence.";
  CharPtr      note = "";

#ifdef WIN_MAC
  bfp = (BaseFormPtr) currentFormDataPtr;
#else
  bfp = (BaseFormPtr) GetObjectExtra (i);
#endif
  update = FALSE;
  if (bfp != NULL) {
    sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
    if (sep != NULL) {
      SeqEntryPack (sep);
      EntryChangeGBSource (sep);
      EntryCheckGBBlock (sep);
      GetRidOfEmptyFeatsDescStrings (0, sep);
      GetRidOfLocusInSeqIds (0, sep);
      if (OkayToWriteTheEntity (bfp->input_entityID, bfp->form, FALSE) == eOkToWriteEntity_Continue) {
        /*SetChecklistValue (checklistForm, 7);*/
        path [0] = '\0';
        StringNCpy_0 (path, bfp->filepath, sizeof (path));
        dfault [0] = '\0';
        GetDefaultTitleFromForm (bfp->form, dfault, sizeof (dfault), bfp->filepath);
        ptr = StringRChr (dfault, '.');
        if (ptr != NULL) {
          *ptr = '\0';
        }
        if (StringLen (dfault) < sizeof (dfault) - 5) {
          StringCat (dfault, ".sqn");
        }
        if (GetOutputFileName (path, sizeof (path), dfault)) {
          update = PropagateFromGenBankBioseqSet (sep, TRUE);
          NormalizeDescriptorOrder (sep);
          update = TRUE; /* because of NormalizeDescriptorOrder */
          if (SeqEntryHasAligns (bfp->input_entityID, sep)) {
            if (Message (MSG_YN, "Remove alignments?") == ANS_YES) {
              SeqEntryExplore (sep, NULL, DoRemoveAlignmentFromRecord);
              update = TRUE;
            }
          }
          SeqMgrClearFeatureIndexes (bfp->input_entityID, NULL);
          if (WriteTheEntityID (bfp->input_entityID, path, FALSE)) {
            email_address = ReturnSubmissionEmailAddress (bfp->input_entityID);
            if (IsAnySequenceMissingAnnotation (sep)) {
              note = missing_annot;
            }
            /*SetChecklistValue (checklistForm, 5);*/
            ptr = StringRChr (path, DIRDELIMCHR);
            if (ptr != NULL) {
              ptr++;
              str = MemNew (sizeof (Char) * (StringLen (ptr) + StringLen (fmt_file) + StringLen (email_address) + StringLen (note)));
              if (str != NULL) {
                sprintf (str, fmt_file, ptr, email_address, note);
                UseWindow ((WindoW) bfp->form);
                Message (MSG_OK, str);
                MemFree (str);
                if (update) {
                  ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
                }
                return;
              }
            }
            str = MemNew (sizeof (Char) * (StringLen (fmt_file) + StringLen (email_address) + StringLen (note)));
            if (str != NULL) {
              sprintf (str, fmt_no_file, email_address, note);
              UseWindow ((WindoW) bfp->form);
              Message (MSG_OK, str);
              MemFree (str);
            }
          } else {
            UseWindow ((WindoW) bfp->form);
            Message (MSG_ERROR, "Unable to write file.");
          }
        } else {
          /*SetChecklistValue (checklistForm, 5);*/
          UseWindow ((WindoW) bfp->form);
        }
      }
    }
    if (update) {
      ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
    }
  }
}

extern Boolean SaveSeqSubmitProc (BaseFormPtr bfp, Boolean saveAs)

{
  Char         dfault [32];
  Char         path [PATH_MAX];
  CharPtr      ptr;
  SeqEntryPtr  sep;
  Char         suffix [32];
  Char         tmp [32];
  Boolean      update;

  if (bfp != NULL) {
    sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
    if (sep != NULL) {
      SeqEntryPack (sep);
      EntryChangeGBSource (sep);
      EntryCheckGBBlock (sep);
      GetRidOfEmptyFeatsDescStrings (0, sep);
      GetRidOfLocusInSeqIds (0, sep);
      path [0] = '\0';
      StringNCpy_0 (path, bfp->filepath, sizeof (path));
      if (StringHasNoText (path) || saveAs) {
        dfault [0] = '\0';
        GetDefaultTitleFromForm (bfp->form, dfault, sizeof (dfault), bfp->filepath);
        ptr = StringRChr (dfault, '.');
        if (ptr != NULL) {
          *ptr = '\0';
        }
        suffix [0] = '\0';
        if (GetSequinAppParam ("PREFERENCES", "SUFFIX", ".sqn", tmp, sizeof (tmp))) {
          if (tmp [0] == '.') {
            StringNCpy_0 (suffix, tmp, sizeof (suffix));
          } else {
            StringCpy (suffix, ".");
            StringNCpy_0 (suffix + 1, tmp, sizeof (suffix) - 1);
          }
        }
        if (StringLen (dfault) < sizeof (dfault) - StringLen (suffix)) {
          StringCat (dfault, suffix);
        }
        if (! (GetOutputFileName (path, sizeof (path), dfault))) return FALSE;
      }
      update = PropagateFromGenBankBioseqSet (sep, TRUE);
      NormalizeDescriptorOrder (sep);
      update = TRUE; /* because of NormalizeDescriptorOrder */
      SeqMgrClearFeatureIndexes (bfp->input_entityID, NULL);
      if (WriteTheEntityID (bfp->input_entityID, path, FALSE)) {
        bfp->filepath = MemFree (bfp->filepath);
        bfp->filepath = StringSave (path);
        if (update) {
          ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
        }
        return TRUE;
      } else {
        Message (MSG_ERROR, "Unable to write file.");
        if (update) {
          ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
        }
        return FALSE;
      }
    }
  }
  return FALSE;
}

static void SaveBinSeqEntry (IteM i)

{
  BaseFormPtr  bfp;
  Char         dfault [32];
  Char         path [PATH_MAX];
  CharPtr      ptr;
  Boolean      saveAs = TRUE;
  SeqEntryPtr  sep;
  Boolean      update;

  bfp = (BaseFormPtr) GetObjectExtra (i);
  if (bfp != NULL) {
    sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
    if (sep != NULL) {
      SeqEntryPack (sep);
      EntryChangeGBSource (sep);
      EntryCheckGBBlock (sep);
      GetRidOfEmptyFeatsDescStrings (0, sep);
      GetRidOfLocusInSeqIds (0, sep);
      path [0] = '\0';
      StringNCpy_0 (path, bfp->filepath, sizeof (path));
      if (StringHasNoText (path) || saveAs) {
        dfault [0] = '\0';
        GetDefaultTitleFromForm (bfp->form, dfault, sizeof (dfault), bfp->filepath);
        ptr = StringRChr (dfault, '.');
        if (ptr != NULL) {
          *ptr = '\0';
        }
        if (StringLen (dfault) < sizeof (dfault) - 5) {
          StringCat (dfault, ".val");
        }
        if (! (GetOutputFileName (path, sizeof (path), dfault))) return;
      }
      update = PropagateFromGenBankBioseqSet (sep, TRUE);
      NormalizeDescriptorOrder (sep);
      update = TRUE; /* because of NormalizeDescriptorOrder */
      SeqMgrClearFeatureIndexes (bfp->input_entityID, NULL);
      if (WriteTheEntityID (bfp->input_entityID, path, TRUE)) {
        bfp->filepath = MemFree (bfp->filepath);
        bfp->filepath = StringSave (path);
        if (update) {
          ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
        }
        return;
      } else {
        Message (MSG_ERROR, "Unable to write file.");
        if (update) {
          ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
        }
      }
    }
  }
}

static CharPtr google_earth_1 =
  "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n" \
  "<kml xmlns=\"http://earth.google.com/kml/2.2\">\n" \
  "<Document>\n" \
  "  <name>KmlFile</name>\n" \
  "  <StyleMap id=\"default_copy0+nicon=http://maps.google.com/mapfiles/kml/pal3/icon60.png+hicon=http://maps.google.com/mapfiles/kml/pal3/icon52.png\">\n" \
  "    <Pair>\n" \
  "      <key>normal</key>\n" \
  "      <styleUrl>#default_copy0+icon=http://maps.google.com/mapfiles/kml/pal3/icon60.png</styleUrl>\n" \
  "    </Pair>\n" \
  "    <Pair>\n" \
  "      <key>highlight</key>\n" \
  "      <styleUrl>#default_copy0+icon=http://maps.google.com/mapfiles/kml/pal3/icon52.png</styleUrl>\n" \
  "    </Pair>\n" \
  "  </StyleMap>\n" \
  "  <Style id=\"default_copy0+icon=http://maps.google.com/mapfiles/kml/pal3/icon60.png\">\n" \
  "    <IconStyle>\n" \
  "      <Icon>\n" \
  "        <href>http://maps.google.com/mapfiles/kml/pal3/icon60.png</href>\n" \
  "      </Icon>\n" \
  "    </IconStyle>\n" \
  "  </Style>\n" \
  "  <Placemark>\n";

static CharPtr google_earth_2 =
  "    <styleUrl>#default_copy0+nicon=http://maps.google.com/mapfiles/kml/pal3/icon60.png+hicon=http://maps.google.com/mapfiles/kml/pal3/icon52.png</styleUrl>\n" \
  "    <Point>\n";

static CharPtr google_earth_3 =
  "    </Point>\n" \
  "  </Placemark>\n" \
  "</Document>\n" \
  "</kml>\n";

static Boolean LaunchedGoogleEarth (
  Uint2 entityID, Uint4 itemID, Uint2 itemtype
)

{
  BioSourcePtr       biop;
  SeqMgrDescContext  context;
  Boolean            format_ok = FALSE;
  FILE               *fp;
  FloatHi            lat = 0.0;
  FloatHi            lon = 0.0;
  CharPtr            lat_lon = NULL;
  Boolean            lat_in_range = FALSE;
  Boolean            lon_in_range = FALSE;
  Char               path [PATH_MAX];
  SeqDescPtr         sdp;
  SubSourcePtr       ssp;
#ifdef OS_UNIX
  Char               cmmd [256];
#endif

  if (itemtype != OBJ_SEQDESC) return FALSE;

  sdp = SeqMgrGetDesiredDescriptor (entityID, NULL, itemID, 0, NULL, &context);
  if (sdp != NULL && sdp->choice == Seq_descr_source) {
    biop = (BioSourcePtr) sdp->data.ptrvalue;
    if (biop != NULL) {
      for (ssp = biop->subtype; ssp != NULL; ssp = ssp->next) {
        if (ssp->subtype != SUBSRC_lat_lon) continue;
        lat_lon = ssp->name;
        if (StringHasNoText (lat_lon)) continue;
        IsCorrectLatLonFormat (lat_lon, &format_ok, &lat_in_range, &lon_in_range);
        if (! format_ok) continue;
        if (! lat_in_range) continue;
        if (! lon_in_range) continue;
        if (! ParseLatLon (lat_lon, &lat, &lon)) continue;
        TmpNam (path);
        /* write to original temp file, so next temp file name will not collide */
        fp = FileOpen (path, "w");
        if (fp != NULL) {
          fprintf (fp, "\n");
          FileClose (fp);
          RememberSqnTempFile (path);
        }
        /* now append .kml extension so proper application is launched */
        StringCat (path, ".kml");
        fp = FileOpen (path, "w");
        if (fp != NULL) {
          fprintf (fp, "%s", google_earth_1);
          fprintf (fp, "    <name>%s</name>\n", lat_lon);
          fprintf (fp, "%s", google_earth_2);
          fprintf (fp, "      <coordinates>%lf,%lf</coordinates>\n", (double) lon, (double) lat);
          fprintf (fp, "%s", google_earth_3);
          FileClose (fp);
          RememberSqnTempFile (path);
#ifdef OS_UNIX
          sprintf (cmmd, "open %s", path);
          system (cmmd);
#endif
#ifdef WIN_MSWIN
          Nlm_MSWin_OpenDocument (path);
#endif
        }
        return TRUE;
      }
    }
  }

  return FALSE;
}

static void LIBCALLBACK ValidNotify (ErrSev sev, int errcode, int subcode,
                                     Uint2 entityID, Uint4 itemID, Uint2 itemtype,
                                     Boolean select, Boolean dblClick, Boolean shftKey)

{
  Int2  handled;

  if (dblClick && entityID > 0 && itemID > 0 && itemtype > 0) {
    if (itemtype == OBJ_SEQDESC && shftKey) {
      if (LaunchedGoogleEarth (entityID, itemID, itemtype)) return;
    }

    WatchCursor ();
    handled = GatherProcLaunch (OMPROC_EDIT, FALSE, entityID, itemID,
                                itemtype, 0, 0, itemtype, 0);
    ArrowCursor ();
    if (handled != OM_MSG_RET_DONE || handled == OM_MSG_RET_NOPROC) {
      Beep ();
    }
  }
}

#ifndef WIN_MAC
static void RemoveUpdateDates (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  ValNodePtr    nextsdp;
  Pointer PNTR  prevsdp;
  ValNodePtr    sdp;

  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sdp = bsp->descr;
    prevsdp = (Pointer PNTR) &(bsp->descr);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sdp = bssp->descr;
    prevsdp = (Pointer PNTR) &(bssp->descr);
  } else return;
  while (sdp != NULL) {
    nextsdp = sdp->next;
    if (sdp->choice == Seq_descr_update_date) {
      *(prevsdp) = sdp->next;
      sdp->next = NULL;
      SeqDescFree (sdp);
    } else {
      prevsdp = (Pointer PNTR) &(sdp->next);
    }
    sdp = nextsdp;
  }
}

static CharPtr stillHasGBMsg =
"Source information in a GenBank block should be in a BioSource.  Should I fixup?";

#ifdef USE_SMARTNET
static void SmartResetProc (IteM i)
{
    BaseFormPtr   bfp;
    ForM          f = NULL;
    Int4          status;
    ObjMgrDataPtr omdp;  
    OMUserDataPtr omudp;
    SMUserDataPtr sm_usr_data;

    if((bfp = (BaseFormPtr) GetObjectExtra (i)) != NULL) {
        f = bfp->form;
        
        omudp = ObjMgrGetUserData(bfp->input_entityID, 0, 0, SMART_KEY);
        omdp = ObjMgrGetData (bfp->input_entityID);
        
        sm_usr_data = (SMUserDataPtr) omudp->userdata.ptrvalue;
        
        if(omdp->dirty == FALSE) {
            status = sm_usr_data->header->status;
            sm_usr_data->header->status = SMStatClosed;
            SMSendMsgToClient(sm_usr_data);
            sm_usr_data->header->status = (SMStatusCode)status;
            return;
        }
    }
    return;
}
#endif

static void LaunchValidatorForDone (BaseFormPtr bfp, SeqEntryPtr sep, FormActnFunc revalProc, FormActnFunc continueProc);
static Boolean SmallInferenceAccnVer (ForM f);
static void ValSeqEntryFormExEx (ForM f, Boolean doAligns, Int2 limit, Boolean inferenceAccnCheck, FormActnFunc revalProc, FormActnFunc doneProc);

#ifdef USE_SMARTNET
static Boolean LIBCALLBACK AllGenBankOrRefSeq (BioseqPtr bsp, SeqMgrBioseqContextPtr bcontext)

{
  SeqMgrDescContext  dcontext;
  MolInfoPtr         mip;
  BoolPtr            resetUpdateDate;
  ValNodePtr         sdp;
  SeqIdPtr           sip;

  resetUpdateDate = (BoolPtr) bcontext->userdata;
  for (sip = bsp->id; sip != NULL; sip = sip->next) {
    if (sip->choice == SEQID_EMBL || sip->choice == SEQID_DDBJ) {
      *resetUpdateDate = FALSE;
    }
  }
  sdp = SeqMgrGetNextDescriptor (bsp, NULL, Seq_descr_molinfo, &dcontext);
  if (sdp != NULL) {
    mip = (MolInfoPtr) sdp->data.ptrvalue;
    if (mip != NULL) {
      if (mip->tech == MI_TECH_htgs_1 || mip->tech == MI_TECH_htgs_2 ||
          mip->tech == MI_TECH_htgs_3 || mip->tech == MI_TECH_htgs_0) {
        *resetUpdateDate = FALSE;
      }
    }
  }
  return TRUE;
}

static void DoSmartReport (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr        bsp;
  BioseqSetPtr     bssp;
  Char             id [42];
  ValNodePtr PNTR  vnpp;

  vnpp = (ValNodePtr PNTR) mydata;
  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    SeqIdWrite (bsp->id, id, PRINTID_FASTA_LONG, sizeof (id) - 1);
    ValNodeCopyStr (vnpp, 0, id);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    ValNodeAddInt (vnpp, bssp->_class, 0);
  }
}

static VoidPtr SmartStructureReport (SeqEntryPtr sep)

{
  ValNodePtr  vnp = NULL;

  SeqEntryExplore (sep, (Pointer) &vnp, DoSmartReport);
  return vnp;
}

static Boolean ValNodeListsDiffer (ValNodePtr vnp1, ValNodePtr vnp2)

{
  while (vnp1 != NULL && vnp2 != NULL) {
    if (vnp1->choice != vnp2->choice) return TRUE;
    if (vnp1->choice == 0) {
      if (StringCmp (vnp1->data.ptrvalue, vnp2->data.ptrvalue) != 0) return TRUE;
    }
    vnp1 = vnp1->next;
    vnp2 = vnp2->next;
  }
  if (vnp1 != NULL || vnp2 != NULL) return TRUE;
  return FALSE;
}

static void SmartnetDoneFuncEx (BaseFormPtr bfp, Boolean validate);
static void SmartnetDoneNoValidateFunc (ForM f);

static void SmartnetDoneValidateFunc (ForM f)
{
  Boolean  inferenceAccnCheck;

  inferenceAccnCheck = SmallInferenceAccnVer (f);
  ValSeqEntryFormExEx (f, TRUE, VALIDATE_ALL, inferenceAccnCheck, SmartnetDoneValidateFunc, SmartnetDoneNoValidateFunc);
}

static void SmartnetDoneNoValidateFunc (ForM f)
{
  BaseFormPtr bfp;

  bfp = (BaseFormPtr) GetObjectExtra (f);
  if (bfp == NULL) return;
  SmartnetDoneFuncEx (bfp, FALSE);
}


static void SmartnetDoneFuncEx (BaseFormPtr bfp, Boolean validate)

{
    MsgAnswer     ans;
    ForM          f = NULL;
    Boolean       hasGBStuff;
    ValNodePtr    sdp;
    SeqEntryPtr   sep;
    Boolean       update;
    
    ObjMgrDataPtr omdp;  
    OMUserDataPtr omudp;
    SMUserDataPtr sm_usr_data;
/*    Uint2         entityID; */

    Boolean       resetUpdateDate = TRUE;
    ValNodePtr    vnp;
    EOkToWriteEntity continue_cancel_validate;

    if(bfp != NULL) {
        f = bfp->form;
        
        omudp = ObjMgrGetUserData(bfp->input_entityID, 0, 0, SMART_KEY);
        if (omudp == NULL) return;

        /* for now set dirty flag to force validation in OkayToWriteTheEntity */
        ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);

        omdp = ObjMgrGetData (bfp->input_entityID);

        sm_usr_data = (SMUserDataPtr) omudp->userdata.ptrvalue;

        if(omdp->dirty == FALSE || (sm_usr_data->header->dirty & 0x02)) {
            sm_usr_data->header->status = SMStatClosed;
            SMSendMsgToClient(sm_usr_data);
            
            /*
            entityID = bfp->input_entityID;
            RemoveSeqEntryViewer (bfp->form);
            ObjMgrFreeUserData(entityID, 0, 0, SMART_KEY); 
            */
            
            /* ObjMgrSendMsg (OM_MSG_DEL, entityID, 0, 0); */
            ObjMgrFree(omdp->datatype, omdp->dataptr);
            
            return;
        }

        sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
        if (sep != NULL) {
            
            if (EntrezASN1Detected (sep)) {
                Message(MSG_OK, 
                        "You may not commit entry retrieved from Entrez.\n"
                        "Please close this window instead");
                return;
            }
            
            SeqEntryPack (sep);
            EntryChangeGBSource (sep);
            hasGBStuff = EntryCheckGBBlock (sep);
            GetRidOfEmptyFeatsDescStrings (0, sep);
            GetRidOfLocusInSeqIds (0, sep);

            if (hasGBStuff) {
                /*  ans = Message (MSG_YNC, stillHasGBMsg);
                    if (ans == ANS_CANCEL) return;
                    if (ans == ANS_YES) { */
                MySeqEntryToAsn3 (sep, TRUE, FALSE, FALSE);
                /* } */
            }
            move_cds (sep);
            /* now instantiating protein titles */
            InstantiateProteinTitles (bfp->input_entityID, NULL);
            update = PropagateFromGenBankBioseqSet (sep, FALSE);
            NormalizeDescriptorOrder (sep);
            update = TRUE; /* because of NormalizeDescriptorOrder */
            SeqMgrClearFeatureIndexes (bfp->input_entityID, NULL);

            if (validate) {
              continue_cancel_validate = OkayToWriteTheEntity (bfp->input_entityID, f, TRUE);
              if (continue_cancel_validate == eOkToWriteEntity_Cancel) {
                  if (update && bfp != NULL) {
                      ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
                  }
                  return;
              } else if (continue_cancel_validate == eOkToWriteEntity_Validate) {
                  /* launch validator */
                  LaunchValidatorForDone (bfp, sep, SmartnetDoneValidateFunc, SmartnetDoneNoValidateFunc);
                  return;
              }
            }
            /* ans = Message (MSG_YN, "Reset Update Date?"); */
            ans = ANS_YES;

            /* reset update date in smart mode if GenBank or RefGene, and not HTGS */

            SeqMgrExploreBioseqs (bfp->input_entityID, 0, (Pointer) &resetUpdateDate, AllGenBankOrRefSeq, TRUE, TRUE, TRUE);
            if (/* ans == ANS_YES */ resetUpdateDate) {
                SeqEntryExplore (sep, NULL, RemoveUpdateDates);
                sdp = CreateNewDescriptor (sep, Seq_descr_update_date);
                if (sdp != NULL) {
                    sdp->data.ptrvalue = DateCurr ();
                }
                PropagateFromGenBankBioseqSet (sep, FALSE);
                NormalizeDescriptorOrder (sep);
            }
            CdCheck (sep, NULL);

            omudp = ObjMgrGetUserData(bfp->input_entityID, 0, 0, DUMB_KEY);
            if (omudp != NULL) {
              vnp = SmartStructureReport (sep);
              if (ValNodeListsDiffer (vnp, omudp->userdata.ptrvalue)) {
                sm_usr_data->header->dirty |= 0x04; /* set rearranged signal */
              }
              ValNodeFreeData (vnp);
            }
            SeqMgrClearFeatureIndexes (bfp->input_entityID, NULL);

        }

        if(sm_usr_data->header->format == OBJ_SEQENTRY) {
            SMWriteBioseqObj(sep, sm_usr_data, NULL);
        } else {
            sm_usr_data->header->format = omdp->datatype;
            SMWriteBioseqObj(omdp->dataptr, sm_usr_data, NULL); 
        }

        /* ObjMgrDelete(omdp->datatype, omdp->dataptr); */

        omdp->dirty = FALSE;
        
        /*
        entityID = bfp->input_entityID;
        RemoveSeqEntryViewer (bfp->form);  
        ObjMgrFreeUserData(entityID, 0, 0, SMART_KEY); 
        */

        HideBioseqView ((WindoW) bfp->form);  

        /* ObjMgrSendMsg (OM_MSG_DEL, entityID, 0, 0); */
        ObjMgrFree(omdp->datatype, omdp->dataptr);

        subtoolRecordDirty = FALSE;
        FileRemove (SEQUIN_EDIT_TEMP_FILE);
        FileRemove (SEQUIN_EDIT_PREV_FILE);
        /* FileRemove (SEQUIN_EDIT_BACK_FILE); */
        FileRemove (SEQUIN_EDIT_ARCH_FILE);
        FileRename (SEQUIN_EDIT_BACK_FILE, SEQUIN_EDIT_ARCH_FILE);
        return;
        
    } else {
        Message(MSG_ERROR, "NULL pointer for brp ...? BUG!BUG!BUG!");
        return;
    }
}


static void SmartnetDoneFunc (BaseFormPtr bfp)
{
  SmartnetDoneFuncEx (bfp, TRUE);
}


static void SmartnetDoneProc (IteM i)

{
	BaseFormPtr bfp;

	bfp = (BaseFormPtr) GetObjectExtra (i);
	SmartnetDoneFunc (bfp);
}
#endif


static void SubtoolDoneFuncEx (ForM f, Boolean validate);
static void SubtoolDoneNoValidateFunc (ForM f);

static void SubtoolDoneValidateFunc (ForM f)
{
  Boolean  inferenceAccnCheck;

  inferenceAccnCheck = SmallInferenceAccnVer (f);
  ValSeqEntryFormExEx (f, TRUE, VALIDATE_ALL, inferenceAccnCheck, SubtoolDoneValidateFunc, SubtoolDoneNoValidateFunc);
}


static void SubtoolDoneNoValidateFunc (ForM f)
{
	SubtoolDoneFuncEx (f, FALSE);
}


static void SubtoolDoneFuncEx (ForM form, Boolean validate)
{
  MsgAnswer    ans;
  BaseFormPtr  bfp;
  ForM         f;
  Boolean      hasGBStuff;
  ValNodePtr   sdp;
  SeqEntryPtr  sep;
  Boolean      update;
  EOkToWriteEntity continue_review_cancel;

  if (subtoolEntityID > 0) {
    sep = GetTopSeqEntryForEntityID (subtoolEntityID);
    if (sep != NULL) {
      SeqEntryPack (sep);
      EntryChangeGBSource (sep);
      hasGBStuff = EntryCheckGBBlock (sep);
      GetRidOfEmptyFeatsDescStrings (0, sep);
      GetRidOfLocusInSeqIds (0, sep);
      if (hasGBStuff) {
        ans = Message (MSG_YNC, stillHasGBMsg);
        if (ans == ANS_CANCEL) return;
        if (ans == ANS_YES) {
          MySeqEntryToAsn3 (sep, TRUE, FALSE, FALSE);
        }
      }
      move_cds (sep);
      /* now instantiating protein titles */
      InstantiateProteinTitles (subtoolEntityID, NULL);
      update = PropagateFromGenBankBioseqSet (sep, FALSE);
      NormalizeDescriptorOrder (sep);
      update = TRUE; /* because of NormalizeDescriptorOrder */
      f = NULL;
      bfp = (BaseFormPtr) GetObjectExtra (form);
      if (bfp != NULL) {
        f = bfp->form;
      }
      SeqMgrClearFeatureIndexes (bfp->input_entityID, NULL);
      if (validate) {
        continue_review_cancel = OkayToWriteTheEntity (subtoolEntityID, f, validate);
        if (continue_review_cancel == eOkToWriteEntity_Cancel) {
          if (update && bfp != NULL) {
            ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
          }
          return;
        } else if (continue_review_cancel == eOkToWriteEntity_Validate) {
          /* launch validator */
          LaunchValidatorForDone (bfp, sep, SubtoolDoneValidateFunc, SubtoolDoneNoValidateFunc);
          return;
        }
      }
      /* ans = Message (MSG_YN, "Reset Update Date?"); */
      ans = ANS_YES;
      if (ans == ANS_YES) {
        SeqEntryExplore (sep, NULL, RemoveUpdateDates);
        sdp = CreateNewDescriptor (sep, Seq_descr_update_date);
        if (sdp != NULL) {
          sdp->data.ptrvalue = DateCurr ();
        }
        PropagateFromGenBankBioseqSet (sep, FALSE);
        NormalizeDescriptorOrder (sep);
      }
      CdCheck (sep, NULL);
    }
    /*SetChecklistValue (checklistForm, 7);*/
    SeqMgrClearFeatureIndexes (bfp->input_entityID, NULL);
    if (WriteTheEntityID (subtoolEntityID, "stdout", FALSE)) {
      subtoolRecordDirty = FALSE;
      FileRemove (SEQUIN_EDIT_TEMP_FILE);
      FileRemove (SEQUIN_EDIT_PREV_FILE);
      /* FileRemove (SEQUIN_EDIT_BACK_FILE); */
      FileRemove (SEQUIN_EDIT_ARCH_FILE);
      FileRename (SEQUIN_EDIT_BACK_FILE, SEQUIN_EDIT_ARCH_FILE);
    } else {
      Message (MSG_POSTERR, "Unable to write ASN.1 file");
      return;
    }
    /*SetChecklistValue (checklistForm, 5);*/
  }
  QuitProgram ();
}


static void SubtoolDoneProc (IteM i)

{
  BaseFormPtr  bfp;

  bfp = (BaseFormPtr) GetObjectExtra (i);
  if (bfp != NULL) {
    SubtoolDoneFuncEx (bfp->form, TRUE);
  }
}
#endif

static Boolean ReviewErrorsForValidationFailure (void)
{
  WindoW w;
  GrouP  h, c;
  PrompT p;
  ButtoN b;
  ModalAcceptCancelData acd;
  
  acd.accepted = FALSE;
  acd.cancelled = FALSE;
  
  w = ModalWindow(-20, -13, -10, -10, NULL);
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);
  
  p = StaticPrompt (h, "Submission has validation errors.", 0, 0, programFont, 'l');
  c = HiddenGroup (h, 3, 0, NULL);
  SetGroupSpacing (c, 10, 10);
  b = PushButton (c, "Review Errors", ModalAcceptButton);
  SetObjectExtra (b, &acd, NULL);
  b = PushButton (c, "Continue", ModalCancelButton);
  SetObjectExtra (b, &acd, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) p, (HANDLE) c, NULL);
  
  Show(w); 
  Select (w);
  while (!acd.accepted && ! acd.cancelled)
  {
    ProcessExternalEvent ();
    Update ();
  }
  ProcessAnEvent ();
  Remove (w);
  if (acd.accepted)
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}


static void LaunchValidatorForDone (BaseFormPtr bfp, SeqEntryPtr sep, FormActnFunc revalProc, FormActnFunc continueProc)
{
  ValidStructPtr  vsp;
  Int2            verbosity = 2;
  Char            buf [32];
  Boolean         allRawOrSeg = TRUE;
  ErrHookProc     oldErrHook;
  ErrSev          oldErrSev;
  Int2            errors;
  Int2            j;

  WatchCursor ();
  Update ();
  vsp = ValidStructNew ();
  if (vsp != NULL) {
    /*SetChecklistValue (checklistForm, 6);*/

    verbosity = 2;
    if (GetSequinAppParam ("SETTINGS", "VALIDATEVERBOSITY", NULL, buf, sizeof (buf))) {
      if (! StrToInt (buf, &verbosity)) {
        verbosity = 2;
      }
    }

    CreateValidateWindowExEx (ValidNotify, "Sequin Validation Errors",
                            programFont, SEV_INFO, verbosity, bfp, revalProc, continueProc, TRUE);
    ClearValidateWindow ();
    SeqEntryExplore (sep, (Pointer) (&allRawOrSeg), CheckForCookedBioseqs);
    if (allRawOrSeg) {
      vsp->useSeqMgrIndexes = TRUE;
    }
    HideValidateDoc ();
    vsp->suppressContext = ShouldSetSuppressContext ();
    vsp->justShowAccession = ShouldSetJustShowAccession ();
    oldErrHook = ErrSetHandler (ValidErrHook);
    oldErrSev = ErrSetMessageLevel (SEV_NONE);
    vsp->validateAlignments = TRUE;
    vsp->alignFindRemoteBsp = TRUE;
    vsp->doSeqHistAssembly = FALSE;
    for (j = 0; j < 6; j++) {
      vsp->errors [j] = 0;
    }
    vsp->errfunc = ValidErrCallback;
    ValidateSeqEntry (sep, vsp);
    if (indexerVersion && useEntrez) {
      TaxonValidate (sep, vsp);
    }
    ErrSetMessageLevel (oldErrSev);
    ErrSetHandler (oldErrHook);
    ErrClear ();
    ShowValidateDoc ();
    errors = 0;
    for (j = 0; j < 6; j++) {
      errors += vsp->errors [j];
    }
    if (errors == 0) {
      ArrowCursor ();
      Message (MSG_OK, "Validation test succeeded.");
      FreeValidateWindow ();
    } else {
      RepopulateValidateFilter ();
    }
    ValidStructFree (vsp);
    /*SetChecklistValue (checklistForm, 5);*/
  }
  ArrowCursor ();
  Update ();
}


static void ProcessDoneButton (ForM f)

{
  Boolean         allRawOrSeg = TRUE;
  BaseFormPtr     bfp;
  Int2            errors;
  Int2            j;
  ErrSev          oldErrSev;
  SeqEntryPtr     sep;
  CharPtr         str;
  ValidStructPtr  vsp;
  CharPtr         fmt_no_file = "Submission is now written.  Please e-mail to %s.%s";
  CharPtr         missing_annot = "  Please include a brief summary of your submission within your correspondence.";
  CharPtr         note = "";
  CharPtr         email_address;

  bfp = (BaseFormPtr) GetObjectExtra (f);
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
#ifndef WIN_MAC
  if (smartnetMode) {
#ifdef USE_SMARTNET
    SmartnetDoneProc ((IteM) f);
#endif
    return;
  }
  if (subtoolMode || stdinMode || binseqentryMode) {
    SubtoolDoneProc ((IteM) f);
    return;
  }
#endif
  WatchCursor ();
  Update ();
  vsp = ValidStructNew ();
  if (vsp != NULL) {
    SeqEntryExplore (sep, (Pointer) (&allRawOrSeg), CheckForCookedBioseqs);
    if (allRawOrSeg) {
      vsp->useSeqMgrIndexes = TRUE;
    }
    if (indexerVersion) {
      vsp->alwaysRequireIsoJTA = TRUE;
    }
    oldErrSev = ErrSetMessageLevel (SEV_MAX);
    vsp->validateAlignments = TRUE;
    vsp->alignFindRemoteBsp = TRUE;
    vsp->doSeqHistAssembly = FALSE;
    for (j = 0; j < 6; j++) {
      vsp->errors [j] = 0;
    }
    ValidateSeqEntry (sep, vsp);
    if (indexerVersion && useEntrez) {
      TaxonValidate (sep, vsp);
    }
    ErrSetMessageLevel (oldErrSev);
    ErrClear ();
    ErrShow ();
    errors = 0;
    for (j = 3; j < 6; j++) {
      errors += vsp->errors [j];
    }
    ValidStructFree (vsp);
    if (errors > 0) {
      ArrowCursor ();
      Update ();
      if (ReviewErrorsForValidationFailure ()) {
        LaunchValidatorForDone (bfp, sep, ProcessDoneButton, NULL);
        return;
      }
    }
    ArrowCursor ();
    Update ();
    if (Message (MSG_YN, "Are you ready to save the record?") == ANS_YES) {
      if (SaveSeqSubmitProc (bfp, TRUE)) {
        if (IsAnySequenceMissingAnnotation (sep)) {
          note = missing_annot;
        }
        email_address = ReturnSubmissionEmailAddress (bfp->input_entityID);
        str = (CharPtr) MemNew (sizeof (Char) * (StringLen (fmt_no_file) + StringLen (email_address) + StringLen (note)));
        if (str != NULL) {
          sprintf (str, fmt_no_file, email_address, note);
          UseWindow ((WindoW) bfp->form);
          Message (MSG_OK, str);
          MemFree (str);
        }
      }
    }
  }
  ArrowCursor ();
  Update ();
}

static void CloseAboutWindowProc (WindoW w)

{
  Remove (w);
}

static void CloseAboutPanelProc (PaneL p, PoinT pt)

{
  WindoW  w;

  w = ParentWindow (p);
  Remove (w);
}

static void AboutProc (IteM i)

{
  PaneL   p;
  WindoW  w;

  w = ModalWindow (-50, -33, -1, -1, CloseAboutWindowProc);
  p = SimplePanel (w, AboutBoxWidth (), AboutBoxHeight (), DrawAbout);
  SetPanelClick (p, NULL, NULL, NULL, CloseAboutPanelProc);
  Show (w);
  Select (w);
}

static void StyleManagerProc (IteM i)

{
  MuskStyleManager ();
}

extern Boolean SequinEntrezInit (CharPtr appl_id, Boolean no_warnings, BoolPtr is_network)

{
  /*
  MonitorPtr  mon;
  Boolean     rsult;

  mon = MonitorStrNewEx ("Sequin", 30, FALSE);
  MonitorStrValue (mon, "Connecting to Entrez service");
  Update ();
  rsult = EntrezInit (appl_id, no_warnings, is_network);
  MonitorFree (mon);
  Update ();
  return rsult;
  */
  return FALSE;
}

/*
#ifndef WIN16
static void Cn3DWinShowProc (IteM i)
{
  WindoW  w;

  if (! BiostrucAvail ()) return;
  if (! EntrezIsInited ()) {
    SequinEntrezInit ("Sequin", FALSE, NULL);
  }
  w = Cn3DWin_Entrez(NULL, useEntrez);
  if (w == NULL) return;
  Show (w);
  Select (w);
}
#endif
*/

typedef struct tax3val {
  Uint2       entityID;
  Uint4       itemID;
  Uint2       itemtype;
  Uint1       organelle;
  OrgRefPtr   orp;
  BioseqPtr   bsp;
  SeqFeatPtr  sfp;
} TaxVal, PNTR TaxValPtr;

typedef struct tax3lst {
  ValNodePtr  head;
  ValNodePtr  tail;
} TaxLst, PNTR TaxLstPtr;

static void RecordSrc (Uint2 entityID, Uint4 itemID, Uint2 itemtype, OrgRefPtr orp,
                       Uint1 organelle, TaxLstPtr tlp, SeqDescrPtr sdp, SeqFeatPtr sfp)

{
  BioseqPtr      bsp;
  BioseqSetPtr   bssp;
  ObjValNodePtr  ovp;
  SeqEntryPtr    sep;
  TaxValPtr      tvp;
  ValNodePtr     vnp;

  if (orp == NULL || tlp == NULL) return;

  tvp = (TaxValPtr) MemNew (sizeof (TaxVal));
  if (tvp == NULL) return;

  vnp = ValNodeNew (tlp->tail);
  if (vnp == NULL) return;

  if (tlp->head == NULL) {
    tlp->head = vnp;
  }
  tlp->tail = vnp;

  tvp->entityID = entityID;
  tvp->itemID = itemID;
  tvp->itemtype = itemtype;
  tvp->organelle = organelle;
  tvp->orp = orp;
  if (sdp != NULL && sdp->extended != 0) {
    ovp = (ObjValNodePtr) sdp;
    if (ovp->idx.parenttype == OBJ_BIOSEQ) {
      bsp = (BioseqPtr) ovp->idx.parentptr;
      if (bsp != NULL) {
        tvp->bsp = bsp;
      }
    } else if (ovp->idx.parenttype == OBJ_BIOSEQSET) {
      bssp = (BioseqSetPtr) ovp->idx.parentptr;
      if (bssp != NULL) {
        sep = bssp->seqentry;
        if (sep != NULL) {
          sep = FindNthBioseq (sep, 1);
          if (sep != NULL) {
            bsp = (BioseqPtr) sep->data.ptrvalue;
            if (bsp != NULL) {
              tvp->bsp = bsp;
            }
          }
        }
      }
    }
  } else if (sfp != NULL) {
    tvp->sfp = sfp;
  }

  vnp->data.ptrvalue = tvp;
}

static void GetSrcDesc (SeqDescrPtr sdp, Pointer userdata)

{
  BioSourcePtr   biop;
  ObjValNodePtr  ovp;
  TaxLstPtr      tlp;

  if (sdp == NULL || sdp->choice != Seq_descr_source) return;
  tlp = (TaxLstPtr) userdata;

  biop = (BioSourcePtr) sdp->data.ptrvalue;
  if (biop == NULL) return;

  if (sdp->extended != 0) {
    ovp = (ObjValNodePtr) sdp;
    RecordSrc (ovp->idx.entityID, ovp->idx.itemID, OBJ_SEQDESC, biop->org, biop->genome, tlp, sdp, NULL);
  }
}

static void GetSrcFeat (SeqFeatPtr sfp, Pointer userdata)

{
  BioSourcePtr  biop;
  TaxLstPtr     tlp;

  if (sfp == NULL || sfp->data.choice != SEQFEAT_BIOSRC) return;
  tlp = (TaxLstPtr) userdata;

  biop = (BioSourcePtr) sfp->data.value.ptrvalue;
  if (biop == NULL) return;

  RecordSrc (sfp->idx.entityID, sfp->idx.itemID, OBJ_SEQFEAT, biop->org, biop->genome, tlp, NULL, sfp);
}

NLM_EXTERN void CDECL  ValidErr VPROTO((ValidStructPtr vsp, int severity, int code1, int code2, const char *fmt, ...));


static void ReportOneBadSpecificHost (ValNodePtr vnp, ValidStructPtr vsp, CharPtr msg_fmt)
{
  ObjValNodePtr ovp;
  BioSourcePtr  biop;
  OrgModPtr     mod;

  if (vnp == NULL || vsp == NULL || StringHasNoText (msg_fmt)) return;

  vsp->sfp = NULL;
  vsp->descr = NULL;
  vsp->bsp = NULL;
  vsp->bssp = NULL;
  biop = NULL;
  mod = NULL;

  if (vnp->choice == OBJ_SEQFEAT)
  {
    vsp->sfp = (SeqFeatPtr) vnp->data.ptrvalue;
    vsp->gcp->entityID = vsp->sfp->idx.entityID;
    vsp->gcp->itemID = vsp->sfp->idx.itemID;
    vsp->gcp->thistype = OBJ_SEQFEAT;
    if (vsp->sfp->idx.parenttype == OBJ_BIOSEQ)
    {
      vsp->bsp = vsp->sfp->idx.parentptr;        
    }
    else if (vsp->sfp->idx.parenttype == OBJ_BIOSEQSET)
    {
      vsp->bssp = vsp->sfp->idx.parentptr;        
    }
    biop = (BioSourcePtr) vsp->sfp->data.value.ptrvalue;
  } 
  else if (vnp->choice == OBJ_SEQDESC)
  {
    vsp->descr = (SeqDescrPtr) vnp->data.ptrvalue;
    if (vsp->descr != NULL && vsp->descr->extended != 0) 
    {
      ovp = (ObjValNodePtr) vsp->descr;
      vsp->gcp->entityID = ovp->idx.entityID;
      vsp->gcp->itemID = ovp->idx.itemID;
      vsp->gcp->thistype = OBJ_SEQDESC;

      if (ovp->idx.parenttype == OBJ_BIOSEQ)
      {
        vsp->bsp = ovp->idx.parentptr;        
      }
      else if (ovp->idx.parenttype == OBJ_BIOSEQSET)
      {
        vsp->bssp = ovp->idx.parentptr;        
      }
    }
    biop = vsp->descr->data.ptrvalue;
  }
  
  if (biop != NULL && biop->org != NULL && biop->org->orgname != NULL)
  {
    mod = biop->org->orgname->mod;
    while (mod != NULL && mod->subtype != ORGMOD_nat_host)
    {
      mod = mod->next;
    }
    if (mod != NULL)
    {      
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_BadSpecificHost, msg_fmt, mod->subname);
    }
  }
}


static void ReportBadSpecificHostValues (SeqEntryPtr sep, ValidStructPtr vsp)
{
  ValNodePtr    misspelled = NULL, bad_caps = NULL, ambiguous = NULL, unrecognized = NULL, vnp;

  Taxon3ValidateSpecificHostsInSeqEntry (sep, &misspelled, &bad_caps, &ambiguous, &unrecognized);

  for (vnp = misspelled; vnp != NULL; vnp = vnp->next) {
    ReportOneBadSpecificHost (vnp, vsp, "Specific host value is misspelled: %s");
  }
  for (vnp = bad_caps; vnp != NULL; vnp = vnp->next) {
    ReportOneBadSpecificHost (vnp, vsp, "Specific host value is incorrectly capitalized: %s");
  } 
  for (vnp = ambiguous; vnp != NULL; vnp = vnp->next) {
    ReportOneBadSpecificHost (vnp, vsp, "Specific host value is ambiguous: %s");
  } 
  for (vnp = unrecognized; vnp != NULL; vnp = vnp->next) {
    ReportOneBadSpecificHost (vnp, vsp, "Invalid value for specific host: %s");
  } 

  misspelled = ValNodeFree (misspelled);
  bad_caps = ValNodeFree (bad_caps);
  unrecognized = ValNodeFree (unrecognized);
}

static Boolean log_tax_asn = FALSE;
static Boolean log_tax_set = FALSE;

static void TaxonValidate (SeqEntryPtr sep, ValidStructPtr vsp)

{
  GatherContext     gc;
  Boolean           has_nucleomorphs;
  Boolean           is_nucleomorph;
  Boolean           is_species_level;
  Boolean           force_tax_consult;
  ValNodePtr        last = NULL;
  OrgRefPtr         orp;
  ErrSev            sev;
  TaxLst            srclist;
  CharPtr           str;
  T3ErrorPtr        t3ep;
  Taxon3RequestPtr  t3rq;
  Taxon3ReplyPtr    t3ry;
  T3DataPtr         tdp;
  T3StatusFlagsPtr  tfp;
  T3ReplyPtr        trp;
  TaxValPtr         tvp;
  ValNodePtr        val;
  ValNodePtr        vnp;
  ValNodePtr        vnp2;

  if (sep == NULL || vsp == NULL) return;
  MemSet ((Pointer) &gc, 0, sizeof (GatherContext));
  vsp->gcp = &gc;

  srclist.head = NULL;
  srclist.tail = NULL;
  VisitDescriptorsInSep (sep, (Pointer) &srclist, GetSrcDesc);
  VisitFeaturesInSep (sep, (Pointer) &srclist, GetSrcFeat);
  if (srclist.head == NULL) return;

  t3rq = Taxon3RequestNew ();
  if (t3rq == NULL) return;

  for (vnp = srclist.head; vnp != NULL; vnp = vnp->next) {
    tvp = (TaxValPtr) vnp->data.ptrvalue;
    if (tvp == NULL) continue;
    orp = AsnIoMemCopy (tvp->orp,
                        (AsnReadFunc) OrgRefAsnRead,
                        (AsnWriteFunc) OrgRefAsnWrite);
    vnp2 = ValNodeAddPointer (&last, 3, (Pointer) orp);
    if (t3rq->request == NULL) {
      t3rq->request = vnp2;
    }
    last = vnp2;
  }

#ifdef OS_UNIX
  if (! log_tax_set) {
    str = (CharPtr) getenv ("LOG_TAX_ASN");
    if (StringDoesHaveText (str)) {
      if (StringICmp (str, "TRUE") == 0) {
        log_tax_asn = TRUE;
      }
    }
    log_tax_set = TRUE;
  }
#endif

  sev = ErrSetMessageLevel (SEV_WARNING);
  if (log_tax_asn) {
    LaunchAsnTextViewer ((Pointer) t3rq, (AsnWriteFunc) Taxon3RequestAsnWrite, "tax3 request");
  }
  t3ry = Tax3SynchronousQuery (t3rq);
  ErrSetMessageLevel (sev);
  Taxon3RequestFree (t3rq);
  if (t3ry == NULL) return;
  if (log_tax_asn) {
    LaunchAsnTextViewer ((Pointer) t3ry, (AsnWriteFunc) Taxon3ReplyAsnWrite, "tax3 result");
  }

  for (trp = t3ry->reply, vnp = srclist.head;
       trp != NULL && vnp != NULL;
       trp = trp->next, vnp = vnp->next) {
    tvp = (TaxValPtr) vnp->data.ptrvalue;
    if (tvp == NULL) continue;
    if (trp->choice == T3Reply_error) {
      t3ep = (T3ErrorPtr) trp->data.ptrvalue;
      if (t3ep != NULL) {
        str = t3ep->message;
        if (str == NULL) {
          str = "?";
        }

        vsp->bssp = NULL;
        vsp->bsp = tvp->bsp;
        vsp->sfp = tvp->sfp;
        vsp->descr = NULL;

        gc.entityID = tvp->entityID;
        gc.itemID = tvp->itemID;
        gc.thistype = tvp->itemtype;

        ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_TaxonomyLookupProblem, "Taxonomy lookup failed with message '%s'", str);
      }
    }
    if (trp->choice != T3Reply_data) continue;
    tdp = (T3DataPtr) trp->data.ptrvalue;
    if (tdp == NULL) continue;

    vsp->bssp = NULL;
    vsp->bsp = tvp->bsp;
    vsp->sfp = tvp->sfp;
    vsp->descr = NULL;

    is_species_level = FALSE;
    has_nucleomorphs = FALSE;
    is_nucleomorph = FALSE;

    for (tfp = tdp->status; tfp != NULL; tfp = tfp->next) {

      /*
      val = tfp->Value_value;
      if (val != NULL && val->choice == Value_value_bool) {
        str = tfp->property;
        if (str == NULL) {
          str = "?";
        }
        if (val->data.intvalue != 0) {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_TaxonomyLookupProblem, "'%s' TRUE", str);
        } else {
          ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_TaxonomyLookupProblem, "'%s' FALSE", str);
        }
      }
      */

      if (StringICmp (tfp->property, "is_species_level") == 0) {
        val = tfp->Value_value;
        if (val != NULL && val->choice == Value_value_bool) {
          is_species_level = (Boolean) (val->data.intvalue != 0);
          if (! is_species_level) {
            gc.entityID = tvp->entityID;
            gc.itemID = tvp->itemID;
            gc.thistype = tvp->itemtype;

            ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_TaxonomyLookupProblem, "Taxonomy lookup reports is_species_level FALSE");
          }
        }
      } else if (StringICmp (tfp->property, "force_consult") == 0) {
        val = tfp->Value_value;
        if (val != NULL && val->choice == Value_value_bool) {
          force_tax_consult = (Boolean) (val->data.intvalue != 0);
          if (force_tax_consult) {
            gc.entityID = tvp->entityID;
            gc.itemID = tvp->itemID;
            gc.thistype = tvp->itemtype;

            ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_TaxonomyLookupProblem, "Taxonomy lookup reports taxonomy consultation needed");
          }
        }
      } else if (StringICmp (tfp->property, "has_nucleomorphs") == 0) {
        val = tfp->Value_value;
        if (val != NULL && val->choice == Value_value_bool) {
          has_nucleomorphs = (Boolean) (val->data.intvalue != 0);
          if (has_nucleomorphs) {
            is_nucleomorph = TRUE;
          }
        }
      }
    }
    if (tvp->organelle == GENOME_nucleomorph && (! is_nucleomorph)) {
      ValidErr (vsp, SEV_WARNING, ERR_SEQ_DESCR_TaxonomyLookupProblem, "Taxonomy lookup does not have expected nucleomorph flag");
    }
  }

  Taxon3ReplyFree (t3ry);
  ValNodeFreeData (srclist.head);

  /* also validate specific-host values */

  ReportBadSpecificHostValues (sep, vsp);
}

static void ValSeqEntryFormExEx (ForM f, Boolean doAligns, Int2 limit, Boolean inferenceAccnCheck, FormActnFunc revalProc, FormActnFunc doneProc)

{
  Boolean         allRawOrSeg = TRUE;
  BaseFormPtr     bfp;
  Char            buf [32];
  Int2            errors;
  Int2            j;
  ErrHookProc     oldErrHook;
  ErrSev          oldErrSev;
  SeqEntryPtr     sep;
  Char            str [32];
  WindoW          validatorWindow = NULL;
  Int2            verbosity;
  ValidStructPtr  vsp;

  bfp = (BaseFormPtr) GetObjectExtra (f);
  if (bfp != NULL) {
    sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
    if (sep != NULL) {
      /*
      if (! EntrezIsInited ()) {
        SequinEntrezInit ("Sequin", FALSE, NULL);
      }
      */
      vsp = ValidStructNew ();
      if (vsp != NULL) {
        WatchCursor ();
        Update ();
        /*SetChecklistValue (checklistForm, 6);*/

       verbosity = 1;
       if (GetSequinAppParam ("SETTINGS", "VALIDATEVERBOSITY", NULL, buf, sizeof (buf))) {
          if (! StrToInt (buf, &verbosity)) {
            verbosity = 1;
          }
        }

        validatorWindow = CreateValidateWindowExExEx (ValidNotify, "Sequin Validation Errors",
                                                  programFont, SEV_INFO, verbosity, bfp,
                                                  revalProc, doneProc, 
                                                  doneProc == NULL && OkToSequester () ? SequesterSequenceList : NULL,
                                                  doneProc == NULL ? SegregateSequenceList : NULL,
                                                  TRUE);
        ClearValidateWindow ();
        SeqEntryExplore (sep, (Pointer) (&allRawOrSeg), CheckForCookedBioseqs);
        if (allRawOrSeg) {
          vsp->useSeqMgrIndexes = TRUE;
        }
        if (indexerVersion) {
          vsp->alwaysRequireIsoJTA = TRUE;
          vsp->farFetchCDSproducts = TRUE;
          vsp->farFetchMRNAproducts = TRUE;
        }
        vsp->validationLimit = limit;
        HideValidateDoc ();
        vsp->suppressContext = ShouldSetSuppressContext ();
        vsp->justShowAccession = ShouldSetJustShowAccession ();
        if (doAligns) {
          vsp->validateAlignments = TRUE;
          vsp->alignFindRemoteBsp = TRUE;
          vsp->doSeqHistAssembly = TRUE;
          vsp->farIDsInAlignments = (Boolean) (subtoolMode || smartnetMode || dirsubMode);
          if (GetSequinAppParam ("SETTINGS", "VALIDATEFARALIGNIDS", NULL, str, sizeof (str))) {
            if (StringICmp (str, "TRUE") == 0) {
              vsp->farIDsInAlignments = TRUE;
            }
          }
        }
        if (useEntrez && inferenceAccnCheck) {
          LookupFarSeqIDs (sep, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE);
          vsp->inferenceAccnCheck = TRUE;
        }
        vsp->testLatLonSubregion = testLatLonSubregion;
        vsp->strictLatLonCountry = strictLatLonCountry;
        vsp->indexerVersion = indexerVersion;
        oldErrHook = ErrSetHandler (ValidErrHook);
        oldErrSev = ErrSetMessageLevel (SEV_NONE);
        for (j = 0; j < 6; j++) {
          vsp->errors [j] = 0;
        }
        vsp->errfunc = ValidErrCallback;
        ValidateSeqEntry (sep, vsp);
        if (indexerVersion && useEntrez) {
          SetTitle (validatorWindow, "Validating Taxonomy");
          Update ();
          TaxonValidate (sep, vsp);
          SetTitle (validatorWindow, "Sequin Validation Errors");
          Update ();
        }
        ErrSetMessageLevel (oldErrSev);
        ErrSetHandler (oldErrHook);
        ErrClear ();
        ShowValidateDoc ();
        errors = 0;
        for (j = 0; j < 6; j++) {
          errors += vsp->errors [j];
        }
        if (errors == 0) {
          ArrowCursor ();
          Message (MSG_OK, "Validation test succeeded.");
          FreeValidateWindow ();
        } else {
          RepopulateValidateFilter ();
        }
        ValidStructFree (vsp);
        /*SetChecklistValue (checklistForm, 5);*/
        ArrowCursor ();
        Update ();
      }
    }
  }
}


static void ValSeqEntryFormEx (ForM f, Boolean doAligns, Int2 limit, Boolean inferenceAccnCheck)

{
  ValSeqEntryFormExEx (f, doAligns, limit, inferenceAccnCheck, ValSeqEntryForm, NULL);
}

static void CountInfAccnVer (SeqFeatPtr sfp, Pointer userdata)

{
  Int4Ptr    countP;
  GBQualPtr  gbq;

  if (sfp == NULL || userdata == NULL) return;
  countP = (Int4Ptr) userdata;

  for (gbq = sfp->qual; gbq != NULL; gbq = gbq->next) {
    if (StringICmp (gbq->qual, "inference") == 0) {
      (*countP)++;
    }
  }
}

static Boolean SmallInferenceAccnVer (ForM f)

{
  BaseFormPtr  bfp;
  Int4         count = 0;
  SeqEntryPtr  sep;

  bfp = (BaseFormPtr) GetObjectExtra (f);
  if (bfp == NULL) return FALSE;

  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return FALSE;

  VisitFeaturesInSep (sep, (Pointer) &count, CountInfAccnVer);

  if (count < 100) return TRUE;

  if (indexerVersion) {
    Message (MSG_POST, "Validation skipping %ld inference accession.version tests",
             (long) count);
  }

  return FALSE;
}

extern void ValSeqEntryForm (ForM f)

{
  Boolean  inferenceAccnCheck;

  inferenceAccnCheck = SmallInferenceAccnVer (f);
  ValSeqEntryFormEx (f, TRUE, VALIDATE_ALL, inferenceAccnCheck);
}

static void ValSeqEntryProc (IteM i)

{
  BaseFormPtr  bfp;
  Boolean      inferenceAccnCheck;

#ifdef WIN_MAC
  bfp = (BaseFormPtr) currentFormDataPtr;
#else
  bfp = (BaseFormPtr) GetObjectExtra (i);
#endif
  if (bfp != NULL) {
    inferenceAccnCheck = SmallInferenceAccnVer (bfp->form);
    ValSeqEntryFormEx (bfp->form, TRUE, VALIDATE_ALL, inferenceAccnCheck);
  }
}

static void ValSeqEntryProcNoAln (IteM i)

{
  BaseFormPtr  bfp;

#ifdef WIN_MAC
  bfp = (BaseFormPtr) currentFormDataPtr;
#else
  bfp = (BaseFormPtr) GetObjectExtra (i);
#endif
  if (bfp != NULL) {
    ValSeqEntryFormEx (bfp->form, FALSE, VALIDATE_ALL, FALSE);
  }
}

static void ValSeqEntryProcInfAccn (IteM i)

{
  BaseFormPtr  bfp;

#ifdef WIN_MAC
  bfp = (BaseFormPtr) currentFormDataPtr;
#else
  bfp = (BaseFormPtr) GetObjectExtra (i);
#endif
  if (bfp != NULL) {
    ValSeqEntryFormEx (bfp->form, FALSE, VALIDATE_ALL, TRUE);
  }
}

static void ValSeqEntryProcSpec (IteM i, Int2 limit)

{
  BaseFormPtr  bfp;

#ifdef WIN_MAC
  bfp = (BaseFormPtr) currentFormDataPtr;
#else
  bfp = (BaseFormPtr) GetObjectExtra (i);
#endif
  if (bfp != NULL) {
    ValSeqEntryFormEx (bfp->form, FALSE, limit, FALSE);
  }
}

static void ValSeqEntryProcInst (IteM i)

{
  ValSeqEntryProcSpec (i, VALIDATE_INST);
}

static void ValSeqEntryProcHist (IteM i)

{
  ValSeqEntryProcSpec (i, VALIDATE_HIST);
}

static void ValSeqEntryProcContext (IteM i)

{
  ValSeqEntryProcSpec (i, VALIDATE_CONTEXT);
}

static void ValSeqEntryProcGraph (IteM i)

{
  ValSeqEntryProcSpec (i, VALIDATE_GRAPH);
}

static void ValSeqEntryProcSet (IteM i)

{
  ValSeqEntryProcSpec (i, VALIDATE_SET);
}

static void ValSeqEntryProcFeat (IteM i)

{
  ValSeqEntryProcSpec (i, VALIDATE_FEAT);
}

static void ValSeqEntryProcDesc (IteM i)

{
  ValSeqEntryProcSpec (i, VALIDATE_DESC);
}

#ifdef USE_SPELL
static void SpellCheckTheForm (ForM f)

{
  Boolean         allRawOrSeg = TRUE;
  BaseFormPtr     bfp;
  Char            buf [32];
  Int2            errors;
  Int2            j;
  MonitorPtr      mon;
  ErrHookProc     oldErrHook;
  ErrSev          oldErrSev;
  SeqEntryPtr     sep;
  Int2            verbosity;
  ValidStructPtr  vsp;

  bfp = (BaseFormPtr) GetObjectExtra (f);
  if (bfp != NULL) {
    sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
    if (sep != NULL) {
      WatchCursor ();
      Update ();
      vsp = ValidStructNew ();
      if (vsp != NULL) {
        /*SetChecklistValue (checklistForm, 6);*/
        mon = MonitorStrNewEx ("SpellCheck", 40, FALSE);
        MonitorStrValue (mon, "Connecting to Spell");
        Update ();
        if (! SpellInit ()) {
          /*SetChecklistValue (checklistForm, 5);*/
          ArrowCursor ();
          MonitorFree (mon);
          Update ();
          Message (MSG_ERROR, "Unable to initialize Spell service.");
          Update ();
          return;
        }
        MonitorStrValue (mon, "Performing Spell Check");
        Update ();
        vsp->spellfunc = SpellCheck;
        vsp->spellcallback = SpellCallBack;
        vsp->onlyspell = TRUE;
        vsp->justwarnonspell = TRUE;

       verbosity = 1;
       if (GetSequinAppParam ("SETTINGS", "VALIDATEVERBOSITY", NULL, buf, sizeof (buf))) {
          if (! StrToInt (buf, &verbosity)) {
            verbosity = 1;
          }
        }

        CreateValidateWindowEx (ValidNotify, "Sequin Spell Check Errors",
                                programFont, SEV_INFO, verbosity, bfp, SpellCheckTheForm, TRUE);
        ClearValidateWindow ();
        SeqEntryExplore (sep, (Pointer) (&allRawOrSeg), CheckForCookedBioseqs);
        if (allRawOrSeg) {
          vsp->useSeqMgrIndexes = TRUE;
        }
        HideValidateDoc ();
        vsp->suppressContext = ShouldSetSuppressContext ();
        vsp->justShowAccession = ShouldSetJustShowAccession ();
        oldErrHook = ErrSetHandler (ValidErrHook);
        oldErrSev = ErrSetMessageLevel (SEV_NONE);
        for (j = 0; j < 6; j++) {
          vsp->errors [j] = 0;
        }
        vsp->errfunc = ValidErrCallback;
        ValidateSeqEntry (sep, vsp);
        ErrSetMessageLevel (oldErrSev);
        ErrSetHandler (oldErrHook);
        ErrClear ();
        ShowValidateDoc ();
        errors = 0;
        for (j = 0; j < 6; j++) {
          errors += vsp->errors [j];
        }
        if (errors == 0) {
          ArrowCursor ();
          Message (MSG_OK, "Spelling check succeeded.");
          FreeValidateWindow ();
        } else {
          RepopulateValidateFilter ();
        }
        ValidStructFree (vsp);
        MonitorStrValue (mon, "Closing Spell Check");
        Update ();
        SpellFini ();
        /*SetChecklistValue (checklistForm, 5);*/
        MonitorFree (mon);
      }
      ArrowCursor ();
      Update ();
    }
  }
}

static void SpellCheckSeqEntryProc (IteM i)

{
  BaseFormPtr  bfp;

#ifdef WIN_MAC
  bfp = (BaseFormPtr) currentFormDataPtr;
#else
  bfp = (BaseFormPtr) GetObjectExtra (i);
#endif
  if (bfp != NULL) {
    SpellCheckTheForm (bfp->form);
  }
}
#endif

extern Int4 MySeqEntryToAsn3Ex (SeqEntryPtr sep, Boolean strip, Boolean correct, Boolean force, Boolean dotaxon);
extern Int4 MySeqEntryToAsn3 (SeqEntryPtr sep, Boolean strip, Boolean correct, Boolean force)

{
  Boolean  dotaxon;

  dotaxon = FALSE;
/*#ifdef INTERNAL_NCBI_SEQUIN*/
  if (indexerVersion) {
    if (subtoolMode || smartnetMode) {
      dotaxon = TRUE;
    }
  }
/*#endif*/
  return MySeqEntryToAsn3Ex (sep, strip, correct, force, dotaxon);
}

static Int2  onthisline;
static Int2  remaining;

static void AddMessageToString (CharPtr str, ValNodePtr vnp)

{
  if (str != NULL && vnp != NULL) {
    StringCat (str, (CharPtr) vnp->data.ptrvalue);
    onthisline++;
    remaining--;
    if (remaining == 0) {
      StringCat (str, "\n\n");
    } else if (onthisline < 3) {
      StringCat (str, ", ");
    } else {
      StringCat (str, ",\n");
      onthisline = 0;
    }
  }
}

static Boolean DisplayTestResults (ValNodePtr head)

{
  MsgAnswer    ans;
  Int2         errors;
  size_t       len;
  CharPtr      str;
  ValNodePtr   vnp;
  Int2         warnings;

  if (head != NULL) {

    vnp = head;
    len = 200;
    while (vnp != NULL) {
      len += StringLen ((CharPtr) vnp->data.ptrvalue) + 2;
      vnp = vnp->next;
    }

    errors = 0;
    warnings = 0;
    vnp = head;
    while (vnp != NULL) {
      if (vnp->choice == 0) {
        warnings++;
      } else {
        errors++;
      }
      vnp = vnp->next;
    }

    str = MemNew (len + 2);
    if (str == NULL) return TRUE;
    str [0] = '\0';

    if (errors > 0) {
      StringCat (str, "The following essential information is missing:\n\n");
      onthisline = 0;
      remaining = errors;
      vnp = head;
      while (vnp != NULL) {
        if (vnp->choice != 0) {
          AddMessageToString (str, vnp);
        }
        vnp = vnp->next;
      }
      StringCat (str, "Please fill in the essential information.");
      ans = Message (MSG_OK, str);
      str = MemFree (str);
      return FALSE;
    }

    if (warnings > 0 && (! indexerVersion)) {
      StringCat (str, "The following desired information is missing:\n\n");
      onthisline = 0;
      remaining = warnings;
      vnp = head;
      while (vnp != NULL) {
        if (vnp->choice == 0) {
          AddMessageToString (str, vnp);
        }
        vnp = vnp->next;
      }
      StringCat (str, "Do you wish to proceed anyway?");
      ans = Message (MSG_YN, str);
      str = MemFree (str);
      if (ans == ANS_NO) return FALSE;
      return TRUE;
    }

    str = MemFree (str);
  }
  return TRUE;
}

extern void JustRegisterSeqEntry (BaseFormPtr bfp, Boolean freeit)

{
  Int2  handled;

  if (bfp != NULL) {
    Hide (bfp->form);
  }
  seqviewprocs.filepath = globalPath;
  seqviewprocs.forceSeparateViewer = TRUE;
  SeqEntrySetScope (NULL);
  handled = GatherProcLaunch (OMPROC_VIEW, FALSE, globalEntityID, 1,
                              OBJ_BIOSEQ, 0, 0, OBJ_BIOSEQ, 0);
  seqviewprocs.filepath = NULL;
  ArrowCursor ();
  if (handled != OM_MSG_RET_DONE || handled == OM_MSG_RET_NOPROC) {
    Message (MSG_FATAL, "Unable to launch viewer.");
  } else {
    SendHelpScrollMessage (helpForm, "Editing the Record", NULL);
  }
  ObjMgrSetOptions (OM_OPT_FREE_IF_NO_VIEW, globalEntityID);
  ObjMgrSetDirtyFlag (globalEntityID, TRUE);
  if (bfp != NULL && freeit) {
    Remove (bfp->form);
  }
}

extern void JustRegisterSeqEntryBtn (ButtoN b)

{
  BaseFormPtr  bfp;

  bfp = (BaseFormPtr) GetObjectExtra (b);
  JustRegisterSeqEntry (bfp, TRUE);
}

static void JustRegisterSeqEntryForm (ForM f)

{
  BaseFormPtr  bfp;

  bfp = (BaseFormPtr) GetObjectExtra (f);
  JustRegisterSeqEntry (bfp, FALSE);
}

extern void AddSubmitBlockToSeqEntry (ForM f)

{
  AffilPtr        affil;
  AuthListPtr     authors;
  BaseFormPtr     bfp;
  CitSubPtr       csp;
  ValNodePtr      head;
  SubmitBlockPtr  sbp;
  SeqSubmitPtr    ssp;
  CitSubPtr       tmp;

  bfp = (BaseFormPtr) GetObjectExtra (f);
  if (bfp != NULL) {
    head = TestForm (bfp->form);
    if (head != NULL) {
      if (! DisplayTestResults (head)) {
        ValNodeFreeData (head);
        return;
      }
      ValNodeFreeData (head);
    }
    Hide (bfp->form);
    /*
    globalsbp = (SequinBlockPtr) FormToPointer (bfp->form);
    if (globalsbp == NULL) {
      Message (MSG_OK, "Record will be a Seq-entry instead of a Seq-submit.");
    }
    */
    sbp = (SubmitBlockPtr) FormToPointer (bfp->form);
    if (sbp == NULL) {
      Message (MSG_OK, "Record will be a Seq-entry instead of a Seq-submit.");
    }
    Update ();
    /*
    globalEntityID = PackageFormResults (globalsbp, globalsep, FALSE);
    */
    globalEntityID = 0;
    if (globalsep != NULL) {
      if (sbp != NULL) {
        if (sbp->contact != NULL && sbp->cit != NULL) {
          tmp = CitSubFromContactInfo (sbp->contact);
          csp = sbp->cit;
          if (csp->authors != NULL) {
            authors = csp->authors;
            if (authors->affil == NULL) {
              if (tmp != NULL && tmp->authors != NULL) {
                authors = tmp->authors;
                affil = authors->affil;
                authors->affil = NULL;
                authors = csp->authors;
                authors->affil = affil;
                if (affil != NULL) {
                  affil->phone = MemFree (affil->phone);
                  affil->fax = MemFree (affil->fax);
                  affil->email = MemFree (affil->email);
                }
              }
            }
          }
          CitSubFree (tmp);
        }
        ssp = SeqSubmitNew ();
        if (ssp != NULL) {
          ssp->datatype = 1;
          ssp->sub = sbp;
          ssp->data = (Pointer) globalsep;
          ObjMgrConnect (OBJ_SEQENTRY, globalsep->data.ptrvalue, OBJ_SEQSUB, (Pointer) ssp);
          if (! ObjMgrRegister (OBJ_SEQSUB, (Pointer) ssp)) {
            ErrPostEx (SEV_ERROR, 0, 0, "ObjMgrRegister failed.");
          }
        } else {
          if (! ObjMgrRegister (OBJ_SEQENTRY, (Pointer) globalsep)) {
            ErrPostEx (SEV_ERROR, 0, 0, "ObjMgrRegister failed.");
          }
        }
      } else {
        if (! ObjMgrRegister (OBJ_SEQENTRY, (Pointer) globalsep)) {
          ErrPostEx (SEV_ERROR, 0, 0, "ObjMgrRegister failed.");
        }
      }
      if (EntrezASN1Detected (globalsep)) {
        ErrPostEx (SEV_WARNING, 0, 0, "This record was retrieved from Entrez.");
      }
      globalEntityID = ObjMgrGetEntityIDForChoice (globalsep);
    }
    globalsbp = NULL;
    globalsep = NULL;
    JustRegisterSeqEntryForm (f);
  }
}

#ifdef WIN_MAC
static void SubmitBlockActivateProc (WindoW w);
static void GenomeFormActivateProc (WindoW w);
#else
#define SubmitBlockActivateProc NULL
#define GenomeFormActivateProc NULL
#endif

CharPtr repackageMsg =
"Do you plan to submit this as an update to one of the databases?";

extern Boolean ProcessOneNucleotideTitle (Int2 seqPackage,
                                          SeqEntryPtr nsep, SeqEntryPtr top);

static void LookForTaxonID (BioSourcePtr biop, Pointer userdata)

{
  DbtagPtr     dbt;
  BoolPtr      notaxid;
  ObjectIdPtr  oip;
  OrgRefPtr    orp;
  ValNodePtr   vnp;

  notaxid = (BoolPtr) userdata;
  if (biop == NULL) return;
  orp = biop->org;
  if (orp == NULL) return;
  for (vnp = orp->db; vnp != NULL; vnp = vnp->next) {
    dbt = (DbtagPtr) vnp->data.ptrvalue;
    if (dbt == NULL) continue;
    if (StringICmp (dbt->db, "taxon") != 0) continue;
    oip = dbt->tag;
    if (oip == NULL) continue;
    if (oip->id != 0) return;
  }
  *notaxid = TRUE;
}

static void RnaProtTrailingCommaFix (SeqFeatPtr sfp, Pointer userdata)

{
  Char        ch;
  size_t      len;
  ProtRefPtr  prp;
  RnaRefPtr   rrp;
  CharPtr     str;
  ValNodePtr  vnp;

  if (sfp == NULL) return;

  if (sfp->data.choice == SEQFEAT_PROT) {
    prp = (ProtRefPtr) sfp->data.value.ptrvalue;
    /* turn trailing space into trailing underscore for validator */
    for (vnp = prp->name; vnp != NULL; vnp = vnp->next) {
      str = (CharPtr) vnp->data.ptrvalue;
      if (StringHasNoText (str)) continue;
      len = StringLen (str);
      if (len < 1) continue;
      ch = str [len - 1];
      while (ch == ' ' && len > 2) {
        len--;
        ch = str [len - 1];
      }
      if (ch == ',') {
        str [len - 1] = '_';
        str [len] = '\0';
      }
    }
  } else if (sfp->data.choice == SEQFEAT_RNA) {
    rrp = (RnaRefPtr) sfp->data.value.ptrvalue;
    /* turn trailing space into trailing underscore for validator */
    if (rrp->ext.choice == 1) {
      str = rrp->ext.value.ptrvalue;
      if (StringDoesHaveText (str)) {
        len = StringLen (str);
        if (len > 0) {
          ch = str [len - 1];
          while (ch == ' ' && len > 2) {
            len--;
            ch = str [len - 1];
          }
          if (ch == ',') {
            str [len - 1] = '_';
            str [len] = '\0';
          }
        }
      }
    }
  }
}

static Boolean HandleOneNewAsnProcEx (BaseFormPtr bfp, Boolean removeold, Boolean askForSubmit,
                                    CharPtr path, Pointer dataptr, Uint2 datatype, Uint2 entityID,
                                    Uint2Ptr updateEntityIDPtr, ValNodePtr PNTR err_list)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  Int2          handled;
  Boolean       notaxid;
  SeqEntryPtr   nsep;
  Boolean       processonenuc;
  SeqEntryPtr   sep;
  ForM          w;

  if (dataptr != NULL && entityID > 0) {
    if (datatype == OBJ_SEQSUB || datatype == OBJ_SEQENTRY ||
        datatype == OBJ_BIOSEQ || datatype == OBJ_BIOSEQSET) {
      WatchCursor ();
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
        VisitFeaturesInSep (sep, NULL, RnaProtTrailingCommaFix);
        /*
        if (seqviewprocs.lockFarComponents) {
          bsplist = LockFarComponents (sep);
        }
        */
        nsep = FindNucSeqEntry (sep);
        processonenuc = TRUE;
        if (IS_Bioseq_set (sep)) {
          bssp = (BioseqSetPtr) sep->data.ptrvalue;
          if (bssp != NULL) {
            if (IsPopPhyEtcSet (bssp->_class)) {
              processonenuc = FALSE;
            } else if (bssp->_class == 7) {
              processonenuc = FALSE;
            } else if (bssp->_class == 1 || bssp->_class == 2) {
              processonenuc = FALSE;
            } else if (bssp->_class == BioseqseqSet_class_gen_prod_set) {
              processonenuc = FALSE;
            } else if (bssp->_class == 0) {
              processonenuc = FALSE;
            } else if (bssp->_class == 255) {
              processonenuc = FALSE;
            }
          }
        }
        if (processonenuc) {
          ProcessOneNucleotideTitle (SEQ_PKG_SINGLE, nsep, sep);
        }
        if (! leaveAsOldAsn) {
          notaxid = FALSE;
          VisitBioSourcesInSep (sep, (Pointer) &notaxid, LookForTaxonID);
          if (notaxid) {
            MySeqEntryToAsn3 (sep, TRUE, FALSE, FALSE);
          }
        }
      }
      if (sep != NULL) {
        if (EntrezASN1Detected (sep)) {
          ErrPostEx (SEV_WARNING, 0, 0, "This record was retrieved from Entrez");
        }
      }
      if ((! indexerVersion) && askForSubmit && datatype != OBJ_SEQSUB) {
        ArrowCursor ();
        Update ();
        if (Message (MSG_YN, repackageMsg) == ANS_YES) {
          /*
          UnlockFarComponents (bsplist);
          */
          globalEntityID = entityID;
          globalsep = sep;
          StringNCpy_0 (globalPath, path, sizeof (globalPath));
          WatchCursor ();
          Update ();
          w = CreateSubmitBlockForm (-50, -33, "Submitting Authors",
                                     FALSE, TRUE, NULL, JustRegisterSeqEntryBtn,
                                     AddSubmitBlockToSeqEntry);
          ArrowCursor ();
          if (w != NULL) {
            Show (w);
            Select (w);
            SendHelpScrollMessage (helpForm, "Submitting Authors Form", NULL);
            return TRUE;
          } else {
            Message (MSG_FATAL, "Unable to create window.");
            return FALSE;
          }
        }
      }
      seqviewprocs.filepath = path;
      seqviewprocs.forceSeparateViewer = TRUE;
      SeqEntrySetScope (NULL);
      handled = GatherProcLaunch (OMPROC_VIEW, FALSE, entityID, 1,
                                  OBJ_BIOSEQ, 0, 0, OBJ_BIOSEQ, 0);
      seqviewprocs.filepath = NULL;
      /*
      UnlockFarComponents (bsplist);
      */
      ArrowCursor ();
      if (handled != OM_MSG_RET_DONE || handled == OM_MSG_RET_NOPROC) {
        Message (MSG_FATAL, "Unable to launch viewer.");
        return FALSE;
      } else {
        SendHelpScrollMessage (helpForm, "Editing the Record", NULL);
      }
      ObjMgrSetOptions (OM_OPT_FREE_IF_NO_VIEW, entityID);
      ObjMgrSetDirtyFlag (entityID, TRUE);
      if (bfp != NULL && removeold) {
        Remove (bfp->form);
      }
      return TRUE;
    } else if (datatype == OBJ_SEQANNOT && dataptr != NULL) {
      entityID = 0;
      if (bfp != NULL) {
        entityID = bfp->input_entityID;
      }
      SeqEntrySetScope (NULL);
      entityID = SmartAttachSeqAnnotToSeqEntry (entityID, (SeqAnnotPtr) dataptr, err_list);
      ArrowCursor ();
      if (entityID != 0) {
        /* code to inhibit multiple updates when attaching multiple feature tables to the same entity */
        if (updateEntityIDPtr != NULL) {
          if (*updateEntityIDPtr == 0) {
            *updateEntityIDPtr = entityID;
            return TRUE;
          } else if (*updateEntityIDPtr == entityID) {
            return TRUE;
          }
        }
        ObjMgrSetDirtyFlag (entityID, TRUE);
        ObjMgrSendMsg (OM_MSG_UPDATE, entityID, 0, 0);
        return TRUE;
      }
    } else if (datatype == OBJ_PROJECT && dataptr != NULL) {
      SeqEntrySetScope (NULL);
      HandleProjectAsn ((ProjectPtr) dataptr, entityID);
      ArrowCursor ();
      return TRUE;
    } else {
      Message (MSG_ERROR, "Unable to process object type %d.", (int) datatype);
      ObjMgrDelete (datatype, dataptr);
    }
  }
  return FALSE;
}

static Boolean HandleOneNewAsnProc (BaseFormPtr bfp, Boolean removeold, Boolean askForSubmit,
                                    CharPtr path, Pointer dataptr, Uint2 datatype, Uint2 entityID,
                                    Uint2Ptr updateEntityIDPtr)
{
  return HandleOneNewAsnProcEx (bfp, removeold, askForSubmit, path, dataptr, datatype, entityID,
                                updateEntityIDPtr, NULL);
}

typedef struct multbioseqform {
  FORM_MESSAGE_BLOCK

  BaseFormPtr  bfp;
  Char         filename [PATH_MAX];
  Boolean      removeold;
  Boolean      askForSubmit;
  ValNodePtr   sephead;
} MultBioseqForm, PNTR MultBioseqFormPtr;

static void CommonHandleMultBioseqs (ButtoN b, Int2 whichbutton)

{
  BioseqPtr          bsp;
  BioseqSetPtr       bssp;
  Uint1              choice = 0;
  Uint2              datatype = 0;
  Uint2              entityID = 0;
  Boolean            is_na = TRUE;
  SeqEntryPtr        list;
  MultBioseqFormPtr  mfp;
  SeqEntryPtr        next;
  ValNodePtr         pip;
  ProjectPtr         proj;
  SeqEntryPtr        sep;

  mfp = (MultBioseqFormPtr) GetObjectExtra (b);
  if (mfp == NULL || mfp->sephead == NULL) return;
  Hide (mfp->form);
  Update ();
  switch (whichbutton) {
    case 1 :
      proj = ProjectNew ();
      if (proj != NULL) {
        pip = ValNodeNew (NULL);
        if (pip != NULL) {
          bsp = (BioseqPtr) mfp->sephead->data.ptrvalue;
          if (bsp != NULL) {
            is_na = (Boolean) ISA_na (bsp->mol);
          }
          if (is_na) {
            choice = ProjectItem_nucent;
          } else {
            choice = ProjectItem_protent;
          }
          pip->choice = choice;
          proj->data = pip;
          pip->data.ptrvalue = (Pointer) mfp->sephead;
          mfp->sephead = NULL;
          HandleProjectAsn (proj, 0);
        }
      }
      break;
    case 2 :
      list = mfp->sephead;
      sep = NULL;
      mfp->sephead = NULL;
      while (list != NULL) {
        next = list->next;
        list->next = NULL;
        if (sep != NULL) {
          AddSeqEntryToSeqEntry (sep, list, TRUE);
        } else {
          sep = list;
        }
        list = next;
      }
      if (sep != NULL) {
        if (sep->choice == 1) {
          datatype = OBJ_BIOSEQ;
        } else if (sep->choice == 2) {
          datatype = OBJ_BIOSEQSET;
        }
        entityID = ObjMgrRegister (datatype, (Pointer) sep->data.ptrvalue);
        HandleOneNewAsnProc (mfp->bfp, mfp->removeold, mfp->askForSubmit,
                             mfp->filename, (Pointer) sep, OBJ_SEQENTRY,
                             entityID, NULL);
      }
      break;
    case 3 :
      sep = ValNodeNew (NULL);
      if (sep != NULL) {
        bssp = BioseqSetNew ();
        if (bssp != NULL) {
          sep->choice = 2;
          sep->data.ptrvalue = (Pointer) bssp;
          bssp->_class = 14;
          bssp->seq_set = mfp->sephead;
          mfp->sephead = NULL;
          SeqMgrSeqEntry (SM_BIOSEQSET, (Pointer) bssp, sep);
          SeqMgrLinkSeqEntry (sep, 0, NULL);
          entityID = ObjMgrRegister (OBJ_SEQENTRY, (Pointer) sep);
          HandleOneNewAsnProc (mfp->bfp, mfp->removeold, mfp->askForSubmit,
                               mfp->filename, (Pointer) sep, OBJ_SEQENTRY,
                               entityID, NULL);
        }
      }
      break;
    default :
      break;
  }
  ArrowCursor ();
  Remove (mfp->form);
  Update ();
}

static void MultToDocSum (ButtoN b)

{
  CommonHandleMultBioseqs (b, 1);
}

static void MultToSegSeq (ButtoN b)

{
  CommonHandleMultBioseqs (b, 2);
}

static void MultToPopSet (ButtoN b)

{
  CommonHandleMultBioseqs (b, 3);
}

static void MultToSingleSeq (ButtoN b)

{
  CommonHandleMultBioseqs (b, 2);
}

static void MultBioseqFormMessage (ForM f, Int2 mssg)

{
  MultBioseqFormPtr  mfp;

  mfp = (MultBioseqFormPtr) GetObjectExtra (f);
  if (mfp != NULL) {
    switch (mssg) {
      case VIB_MSG_CLOSE :
        Remove (f);
        break;
      default :
        if (mfp->appmessage != NULL) {
          mfp->appmessage (f, mssg);
        }
        break;
    }
  }
}

static Boolean ProcessMultipleBioseqs (BaseFormPtr bfp, CharPtr filename, Boolean removeold,
                                    Boolean askForSubmit, ValNodePtr sephead)

{
  ButtoN             b;
  GrouP              c;
  MultBioseqFormPtr  mfp;
  StdEditorProcsPtr  sepp;
  WindoW             w;
#ifndef WIN_MAC
  MenU               m;
#endif

  mfp = (MultBioseqFormPtr) MemNew (sizeof (MultBioseqForm));
  if (mfp == NULL) return FALSE;

  if (!FixIDsAndTitles (sephead, NULL, TRUE)) {
    return FALSE;
  }

  w = FixedWindow (-50, -33, -10, -10, "Sequence Input", NULL);
  SetObjectExtra (w, mfp, StdCleanupFormProc);
  mfp->form = (ForM) w;
  mfp->formmessage = MultBioseqFormMessage;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    SetActivate (w, sepp->activateForm);
    mfp->appmessage = sepp->handleMessages;
  }

  mfp->bfp = bfp;
  mfp->removeold = removeold;
  mfp->askForSubmit = askForSubmit;
  mfp->sephead = sephead;
  StringNCpy_0 (mfp->filename, filename, sizeof (mfp->filename));

#ifndef WIN_MAC
  m = PulldownMenu (w, "File");
  FormCommandItem (m, "Close", (BaseFormPtr) mfp, VIB_MSG_CLOSE);
#endif

  c = HiddenGroup (w, 0, 3, NULL);
  SetGroupSpacing (c, 10, 10);
  b = PushButton (c, "Load into Document Window", MultToDocSum);
  SetObjectExtra (b, mfp, NULL);
  if (sephead->next != NULL) {
    b = PushButton (c, "Load as Segmented Sequence", MultToSegSeq);
    SetObjectExtra (b, mfp, NULL);
    b = PushButton (c, "Load as Population Study", MultToPopSet);
    SetObjectExtra (b, mfp, NULL);
  } else {
    b = PushButton (c, "Load as Single Sequence", MultToSingleSeq);
    SetObjectExtra (b, mfp, NULL);
  }

  RealizeWindow (w);
  Show (w);
  Select (w);
  ArrowCursor ();
  Update ();
  return TRUE;
}

static void ProcessMultipleSimpleSeqs (BaseFormPtr bfp, CharPtr filename, Boolean removeold,
                                       Boolean askForSubmit, ValNodePtr simples)

{
  EntrezGlobalsPtr  egp;

  egp = (EntrezGlobalsPtr) GetAppProperty ("EntrezGlobals");
  if (egp == NULL || egp->retrieveSimpleProc == NULL) return;
  egp->retrieveSimpleProc (NULL, simples);
}

static Boolean HandledAnnotatedProteins (BaseFormPtr bfp, ValNodePtr bioseqs)

{
  BioseqPtr    bsp;
  Int2         code;
  ValNodePtr   descr;
  MolInfoPtr   mip;
  BioseqPtr    nucbsp;
  ValNodePtr   sdp;
  SeqEntryPtr  sep;
  SeqLocPtr    slp;
  CharPtr      title;
  SeqEntryPtr  top = NULL;
  ValNode      vn;
  ValNodePtr   vnp;

  if (bfp == NULL || bioseqs == NULL) return FALSE;
  nucbsp = GetBioseqGivenIDs (bfp->input_entityID, bfp->input_itemID, bfp->input_itemtype);
  if (nucbsp == NULL) return FALSE;
  if (! ISA_na (nucbsp->mol)) return FALSE;
  /* top = GetBestTopParentForData (bfp->input_entityID, nucbsp); */
  top = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (top == NULL) return FALSE;
  for (vnp = bioseqs; vnp != NULL; vnp = vnp->next) {
    bsp = (BioseqPtr) vnp->data.ptrvalue;
    if (bsp == NULL) return FALSE;
    if (! ISA_aa (bsp->mol)) return FALSE;
    title = BioseqGetTitle (bsp);
    if (title == NULL) return FALSE;
    if (StringChr (title, '[') == NULL) return FALSE;
  }
  code = SeqEntryToGeneticCode (top, NULL, NULL, 0);
  SetBatchSuggestNucleotide (nucbsp, code);
  descr = ExtractBioSourceAndPubs (top);
  vn.choice = SEQLOC_WHOLE;
  vn.data.ptrvalue = (Pointer) nucbsp->id;
  slp = &vn;
  for (vnp = bioseqs; vnp != NULL; vnp = vnp->next) {
    bsp = (BioseqPtr) vnp->data.ptrvalue;
    bsp->id = SeqIdFree (bsp->id);
    bsp->id = MakeNewProteinSeqId (slp, NULL);
    SeqMgrReplaceInBioseqIndex (bsp);
    sep = SeqMgrGetSeqEntryForData (bsp);
    if (sep != NULL) {
      mip = MolInfoNew ();
      if (mip != NULL) {
        mip->biomol = 8;
        mip->tech = 13;
        sdp = CreateNewDescriptor (sep, Seq_descr_molinfo);
        if (sdp != NULL) {
          sdp->data.ptrvalue = (Pointer) mip;
        }
      }
      AddSeqEntryToSeqEntry (top, sep, TRUE);
      AutomaticProteinProcess (top, sep, code, FALSE, NULL);
      ValNodeExtract (&(bsp->descr), Seq_descr_title);
    }
  }
  ClearBatchSuggestNucleotide ();
  ReplaceBioSourceAndPubs (top, descr);
  return TRUE;
}

static Boolean DoReadAnythingLoop (BaseFormPtr bfp, CharPtr filename, CharPtr path,
                                   Boolean removeold, Boolean askForSubmit, Boolean alwaysMult,
                                   Boolean parseFastaSeqId, Boolean fastaAsSimpleSeq)

{
  ValNodePtr     bioseqs;
  BioseqPtr      bsp;
  Pointer        dataptr;
  Uint2          datatype;
  Boolean        each;
  Uint2          entityID;
  FILE           *fp;
  ValNodePtr     head = NULL;
  ValNodePtr     last = NULL;
  OMUserDataPtr  omudp;
  ValNodePtr     projects;
  Boolean        rsult;
  SeqEntryPtr    sep;
  SeqEntryPtr    sephead = NULL;
  ValNodePtr     simples;
  Uint2          updateEntityID;
  ValNodePtr     vnp, err_list, vnp_err;
  LogInfoPtr     lip;

  if (filename == NULL) return FALSE;
  fp = FileOpen (filename, "r");
  if (fp != NULL) {
    rsult = FALSE;
    while ((dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, NULL, FALSE, FALSE,
                                              parseFastaSeqId, fastaAsSimpleSeq)) != NULL) {
      vnp = ValNodeAddPointer (&last, datatype, dataptr);
      if (head == NULL) {
        head = vnp;
      }
      last = vnp;
    }
    FileClose (fp);
    bioseqs = ValNodeExtractList (&head, OBJ_BIOSEQ);
    projects = ValNodeExtractList (&head, OBJ_PROJECT);
    simples = ValNodeExtractList (&head, OBJ_FASTA);
    updateEntityID = 0;
    lip = OpenLog ("Errors");
    for (vnp = head; vnp != NULL; vnp = vnp->next) {
      datatype = vnp->choice;
      dataptr = vnp->data.ptrvalue;
      entityID = ObjMgrRegister (datatype, dataptr);
      
      err_list = NULL;
      each = HandleOneNewAsnProcEx (bfp, removeold, askForSubmit, path, dataptr, datatype, entityID, &updateEntityID, &err_list);
      if (err_list != NULL) {
        for (vnp_err = err_list; vnp_err != NULL; vnp_err = vnp_err->next) {
          fprintf (lip->fp, "%s\n", (CharPtr) vnp_err->data.ptrvalue);
          lip->data_in_log = TRUE;
        }
        err_list = ValNodeFreeData (err_list);
        if (datatype == OBJ_SEQANNOT) {
          ExportSeqAnnotFeatureTable (lip->fp, dataptr);
        }
      }
      removeold = FALSE;
      rsult = (Boolean) (rsult || each);
      if (backupMode) {
        subtoolEntityID = entityID;
        omudp = ObjMgrAddUserData (subtoolEntityID, 0, 0, 0);
        if (omudp != NULL) {
          omudp->messagefunc = BackupModeMsgFunc;
        }
        subtoolRecordDirty = TRUE;
      }      
    }
    CloseLog (lip);
    lip = FreeLog (lip);
    ValNodeFree (head);
    if (updateEntityID != 0) {
      ObjMgrSetDirtyFlag (updateEntityID, TRUE);
      ObjMgrSendMsg (OM_MSG_UPDATE, updateEntityID, 0, 0);
    }
    for (vnp = projects; vnp != NULL; vnp = vnp->next) {
      HandleProjectAsn ((ProjectPtr) vnp->data.ptrvalue, 0);
    }
    ValNodeFree (projects);
    if (bioseqs != NULL) {
      if (HandledAnnotatedProteins (bfp, bioseqs)) {
        rsult = TRUE;
        ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
        ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
      } else if (bioseqs->next == NULL && (! alwaysMult)) {
        bsp = (BioseqPtr) bioseqs->data.ptrvalue;
        bioseqs->data.ptrvalue = NULL;
        if (bsp != NULL) {
          sep = SeqMgrGetSeqEntryForData (bsp);
          if (sep == NULL) {
            sep = SeqEntryNew ();
            if (sep != NULL) {
              sep->choice = 1;
              sep->data.ptrvalue = bsp;
              SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) bsp, sep);
            }
          }
        }
        entityID = ObjMgrRegister (OBJ_BIOSEQ, (Pointer) bsp);
        each = HandleOneNewAsnProc (bfp, FALSE, FALSE, path, (Pointer) bsp, OBJ_BIOSEQ, entityID, NULL);
        removeold = FALSE;
        rsult = (Boolean) (rsult || each);
      } else {
        sephead = NULL;
        for (vnp = bioseqs; vnp != NULL; vnp = vnp->next) {
          bsp = (BioseqPtr) vnp->data.ptrvalue;
          if (bsp != NULL) {
            sep = SeqMgrGetSeqEntryForData (bsp);
            if (sep == NULL) {
              sep = SeqEntryNew ();
              if (sep != NULL) {
                sep->choice = 1;
                sep->data.ptrvalue = bsp;
                SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) bsp, sep);
              }
            }
            if (sep != NULL) {
              ValNodeLink (&sephead, sep);
            }
          }
        }
        if (ProcessMultipleBioseqs (bfp, filename, removeold, askForSubmit, sephead)) {
          bioseqs = NULL;
          rsult = TRUE;
        } else {
          for (vnp = bioseqs; vnp != NULL; vnp = vnp->next) {
            bsp = (BioseqPtr) vnp->data.ptrvalue;
            if (bsp != NULL) {
              bsp->idx.deleteme = TRUE;
            }
          }
          DeleteMarkedObjects (0, OBJ_SEQENTRY, sephead);
        }
      }
    }
    ValNodeFree (bioseqs);
    if (simples != NULL) {
      ProcessMultipleSimpleSeqs (bfp, filename, removeold, askForSubmit, simples);
      rsult = TRUE;
    }
    ValNodeFree (simples);
    return rsult;
  }
  return FALSE;
}

/* Sarah's stuff */

#define SQACT_TXTALN 1
#define SQACT_NOTTXTALN 2
#define SQACT_UNPRINTABLECHARS 3

static Int4 SQACT_GuessWhatIAm (FILE *fp)
{
   Boolean  found;
   Int4     len;
   CharPtr  line;
   Int4     prevlen;
   Int4     seq;
   CharPtr  tmp;
   ValNodePtr current_data = NULL;
   CharPtr    cp;

   if (fp == NULL) return SQACT_NOTTXTALN;
   line = MyFGetLine (fp, &current_data);
   found = FALSE;
   while (line != NULL && !found) /* find first non-empty line */
   {
      if (!StringHasNoText(line))
      {
         if (StringLen(line) > 0)
            found = TRUE;
      } 
      else
      {
         MemFree(line);
         line = MyFGetLine (fp, &current_data);
      }
   }
   if (!found) 
   {
     FreeBufferedReadList (current_data);
     return -1;
   }

   /* look for unprintable characters */      
   for (cp = line; *cp != 0; cp++)
   {
      if ((Int4)(*cp) < 0 || (Int4)(*cp) >= 256 || ! isprint ((Int4)(*cp)) && ! isspace ((Int4)(*cp)))
      {
         FreeBufferedReadList (current_data);
         return SQACT_UNPRINTABLECHARS;
      }
   }

   tmp = StringStr(line, "::="); /* ASN.1 */
   if (tmp != NULL)
   {
      MemFree(line);
      line = NULL;
      FreeBufferedReadList (current_data);
      return SQACT_NOTTXTALN;
   }
   tmp = StringStr(line, "LOCUS");  /* GenBank flat file */
   if (tmp != NULL)
   {
      MemFree(line);
      line = NULL;
      FreeBufferedReadList (current_data);
      return SQACT_NOTTXTALN;
   }
   if (line[0] == '>')  /* FASTA sequences or FASTA alignment */
   {
      if (StringNCmp (line, ">PubMed", 7) == 0 ||
          StringNCmp (line, ">Protein", 8) == 0 ||
          StringNCmp (line, ">Nucleotide", 11) == 0 ||
          StringNCmp (line, ">Structure", 10) == 0 ||
          StringNCmp (line, ">Genome", 7) == 0 ||
          StringNCmp (line, ">Feature", 8) == 0 ||
          StringNCmp (line, ">Vector", 7) == 0 ||
          StringNCmp (line, ">Restriction", 12) == 0 ||
          StringNCmp (line, ">Contig", 7) == 0 ||
          StringNCmp (line, ">Virtual", 8) == 0 ||
          StringNCmp (line, ">Message", 8) == 0) {
        MemFree(line);
        line = NULL;
        FreeBufferedReadList (current_data);
      	return SQACT_NOTTXTALN;
      }
      MemFree(line);
      line = MyFGetLine (fp, &current_data);
      prevlen = -1;
      len = 0;
      seq = 0;
      while (line != NULL)
      {
         seq++;
         while (line != NULL && line[0] != '>')
         {
            tmp = StringStr(line, "-");
            if (tmp != NULL)  /* found a gap -> alignment */
            {
               MemFree(line);
               line = NULL;
               FreeBufferedReadList (current_data);
               return SQACT_TXTALN;
            }
            len += StringLen(line);
            MemFree(line);
            line = MyFGetLine (fp, &current_data);
         }
         if (line != NULL)
         {
            MemFree(line);
            line = MyFGetLine (fp, &current_data);
         }
         if (prevlen == -1)
            prevlen = len;
         if (len != prevlen)  /* sequences of different length -> FASTA seq */
         {
            MemFree(line);
            line = NULL;
            FreeBufferedReadList (current_data);
            return SQACT_NOTTXTALN;
         }
      }
      MemFree (line);
      line = NULL;
      FreeBufferedReadList (current_data);
      current_data = NULL;
      if (seq > 1) 
         return SQACT_TXTALN;
      else
         return SQACT_NOTTXTALN;
   } else
   {
      MemFree(line);
      line = NULL;
      FreeBufferedReadList (current_data);
      return SQACT_TXTALN;  /* NEXUS, PHYLIP, etc */
   }
}

static Uint1 sq_transform_code(Uint1 res, Int4 which)
{
   if (which == 4) /* ncbi2na */
   {
      if (res == 1)
         return 0;
      else if (res == 3)
         return 1;
      else if (res == 7)
         return 2;
      else if (res == 18)
         return 3;
      else
         return 0;
   } else if (which == 2) /* ncbi4na */
   {
      if (res == 1)
         return 1;
      else if (res == 3)
         return 2;
      else if (res == 12)
         return 3;
      else if (res == 7)
         return 4;
      else if (res == 16)
         return 5;
      else if (res == 17)
         return 6;
      else if (res == 19)
         return 7;
      else if (res == 18)
         return 8; 
      else if (res == 20)
         return 9;
      else if (res == 22)
         return 10;
      else if (res == 8)
         return 11;
      else if (res == 10)
         return 12;
      else if (res == 4)
         return 13;
      else if (res == 2)
         return 14;
      else if (res == 13)
         return 15;
       else
         return 15;
   }
   return 0;
}

static void SQACT_FixBioseqs (BioseqPtr bsp, Pointer userdata)
{
   Uint1Ptr           array;
   ByteStorePtr       bs;
   Int4               compact;
   SeqMgrDescContext  context;
   Int4               i;
   Int4               j;
   Int4               len;
   MolInfoPtr         mip;
   Uint1              res;
   Uint1              res1;
   Int4               shift;
   SeqPortPtr         spp;
   ValNodePtr         vnp;

   if (bsp == NULL)
      return;
   if (bsp->seq_data_type == Seq_code_gap) return;
   len = bsp->length-1;
   compact = 2;
   spp = SeqPortNew(bsp, 0, len, 0, Seq_code_ncbistdaa);
   while ((res = SeqPortGetResidue(spp)) != SEQPORT_EOF)
   {
      if (res != 1 && res != 3 && res != 0 && res != 18 && res != 7 && res != 13 && res != 24 && res != 21)
      {
         SeqPortFree(spp);
         return;
      } else if (res != 1 && res != 3 && res != 7 && res != 18)
         compact = 4;
   }
   SeqPortFree(spp);
   bs = BSNew((bsp->length+compact-1)/compact);
   array = (Uint1Ptr)MemNew((bsp->length+compact)*sizeof(Uint1));
   spp = SeqPortNew(bsp, 0, len, 0, Seq_code_ncbistdaa);
   for (i=0; i<bsp->length; i++)
   {
      res = SeqPortGetResidue(spp);
      array[i] = sq_transform_code(res, compact);
   }
   for (i=0; i<compact; i++)
      array[i+bsp->length] = 0;
   if (compact == 4)
      shift = 2;
   else
      shift = 4;
   for (j=0; j<bsp->length; j+=compact)
   {
      res = 0;
      res1 = 0;
      for (i=0; i<compact; i++)
      {
         res = array[j+i];
         res1<<=shift;
         res1 |= res;
      }
      BSPutByte(bs, res1);
   }
   MemFree(bsp->seq_data);
   bsp->seq_data = (SeqDataPtr) bs;
   if (compact == 4)
      bsp->seq_data_type = Seq_code_ncbi2na;
   else
      bsp->seq_data_type = Seq_code_ncbi4na;
   bsp->mol = Seq_mol_na;
   bsp->repr = Seq_repr_raw;
   SeqMgrIndexFeatures(0, (Pointer)(bsp));
   vnp = NULL;
   vnp = SeqMgrGetNextDescriptor(bsp, vnp, Seq_descr_molinfo, &context);
   if (vnp != NULL)
   {
      vnp = context.sdp;
      mip = (MolInfoPtr)(vnp->data.ptrvalue);
      mip->biomol = 1;
   }
   MemFree(array);
}

static Boolean CommonReadNewAsnProc (Handle obj, Boolean removeold, Boolean askForSubmit)

{
  AsnIoPtr     aip;
  BaseFormPtr  bfp;
  Pointer      dataptr;
  Uint2        datatype;
  Uint2        entityID;
  FILE         *fp;
  Char         path [PATH_MAX];
  Boolean      rsult;
  SeqEntryPtr  sep;
  Int4         type;
  Uint2        updateEntityID;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (obj);
#endif
  rsult = FALSE;
  path [0] = '\0';
  if (path [0] != '\0' || GetInputFileName (path, sizeof (path), "", "TEXT")) {
    Update ();
    fp = FileOpen (path, "r");
    if (fp == NULL) {
      ArrowCursor ();
      ErrPostEx (SEV_WARNING, 0, 0, "Unable to open file '%s'", path);
      return FALSE;
    }
    type = SQACT_GuessWhatIAm (fp);
    FileClose (fp);
    if (type == SQACT_UNPRINTABLECHARS || type == SQACT_TXTALN)
    {
      sep = NULL;
      if (type == SQACT_TXTALN)
      {
        sep = ReadAnyAlignment (TRUE, path);
      }
      if (sep == NULL)
      {
        aip = AsnIoOpen (path, "rb");
        if (aip != NULL) {
          sep = SeqEntryAsnRead (aip, NULL);
          AsnIoClose (aip);
        }
        if (sep == NULL)
        {
          ArrowCursor ();
          if (type == SQACT_UNPRINTABLECHARS)
          {
            ErrPostEx (SEV_WARNING, 0, 0, "File '%s' contains unprintable characters.  If this is a Word document, you need to save it as plain text.", path);
          }
          else
          {
            ErrPostEx (SEV_WARNING, 0, 0, "Error reading file");
          }
          return FALSE;
        }        
      }
      VisitBioseqsInSep (sep, NULL, SQACT_FixBioseqs);
      if (IS_Bioseq (sep)) {
        datatype = OBJ_BIOSEQ;
      } else {
        datatype = OBJ_BIOSEQSET;
      }
      dataptr = (Pointer) sep->data.ptrvalue;
      entityID = ObjMgrRegister (datatype, dataptr);
      rsult = HandleOneNewAsnProc (bfp, FALSE, askForSubmit, path, dataptr, datatype, entityID, &updateEntityID);
      if (updateEntityID != 0) 
      {
        ObjMgrSetDirtyFlag (updateEntityID, TRUE);
        ObjMgrSendMsg (OM_MSG_UPDATE, updateEntityID, 0, 0);
      }
    } else {
      rsult = DoReadAnythingLoop (bfp, path, path, removeold, askForSubmit, FALSE, TRUE, FALSE);
    }
  }
  ArrowCursor ();
  return rsult;
}

static void FinishAlignmentRead (Handle obj, SeqEntryPtr sep, CharPtr path) {
  BaseFormPtr  bfp;
  Pointer      dataptr;
  Uint2        datatype;
  Uint2        entityID;
  Boolean      rsult;
  Uint2        updateEntityID;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (obj);
#endif

  if (sep == NULL)
  {
    return;
  }
  VisitBioseqsInSep (sep, NULL, SQACT_FixBioseqs);
  if (IS_Bioseq (sep)) {
    datatype = OBJ_BIOSEQ;
  } else {
    datatype = OBJ_BIOSEQSET;
  }
  dataptr = (Pointer) sep->data.ptrvalue;
  entityID = ObjMgrRegister (datatype, dataptr);
  rsult = HandleOneNewAsnProc (bfp, FALSE, FALSE, path,
                               dataptr, datatype, entityID, &updateEntityID);
  if (updateEntityID != 0) {
    ObjMgrSetDirtyFlag (updateEntityID, TRUE);
    ObjMgrSendMsg (OM_MSG_UPDATE, updateEntityID, 0, 0);
  }
}


static void PrintExtraErrorInstructions (FILE *fp, CharPtr message)
{
  CharPtr explanation, end;
  Char    tmp = '\0';
  if (fp == NULL || message == NULL) return;

  if (StringStr (message, "bad characters") == NULL
      && StringStr (message, " found at position ") == NULL) {
    return;
  }

  explanation = StringRChr (message, '(');
  if (explanation == NULL) return;
  if (StringNCmp (explanation, "(expect only ", 13) == 0) {
    end = StringStr (explanation + 13, " here)");
    if (end != NULL) {
      tmp = *end;
      *end = 0;
    }
    fprintf (fp, 
             "Try changing the sequence character specifications for %s.\n",
             explanation + 13);
    if (StringNCmp (explanation + 13, "beginning", 9) == 0) {
      fprintf (fp, 
"\nWhen some of the sequences in an alignment are shorter or "
"longer than others, beginning gap characters are added to "
"the beginning of the sequence to maintain the correct spacing."
"  These will not appear in your sequence file.\n");
    } else if (StringNCmp (explanation + 13, "end", 3) == 0) {
      fprintf (fp,
"\nWhen some of the sequences in an alignment are shorter or "
"longer than others, end gap characters are added to "
"the end of the sequence to maintain the correct spacing."
"  These will not appear in your sequence file.\n");
    } else {
      fprintf (fp,
"\nMiddle gap characters are used to maintain the spacing "
"inside an alignment.  These are not nucleotides and will "
"not appear as part of your sequence file.\n"
"Ambiguous/unknown characters are used to represent indeterminate/ambiguous "
"nucleotides.  These will appear in your sequence file as 'n'.\n"
"Match characters are used to indicate positions where "
"sequences are identical to the first sequence.  These will be "
"replaced by the actual characters from the first sequence.\n");
    }
    if (end != NULL) {
      *end = tmp;
    }
  } else if (StringCmp (explanation, 
                        "(can't specify match chars in first sequence).") == 0) {
    fprintf (fp, "Try changing the match character specification.\n");
    fprintf (fp,
"\nMatch characters are used to indicate positions where "
"sequences are identical to the first sequence.  These will be "
"replaced by the actual characters from the first sequence.\n");
  }
}
static void PrintError (FILE *fp, TErrorInfoPtr eip)
{
  if (eip == NULL || fp == NULL) return;

  fprintf (fp, "*****\nError category %d\n", eip->category);
  if (eip->line_num > -1) {
    fprintf (fp, "Line number %d\n", eip->line_num);
  }
  if (eip->id != NULL) {
    fprintf (fp, "Sequence ID %s\n", eip->id);
  }
  if (eip->message != NULL) {
    fprintf (fp, "%s\n", eip->message);
    PrintExtraErrorInstructions (fp, eip->message);
  }
}
  
static void WalkErrorList (TErrorInfoPtr list, FILE *fp)
{
  TErrorInfoPtr eip;
  
  if (list == NULL || fp == NULL) return;

  for (eip = list; eip != NULL; eip = eip->next) {
    PrintError (fp, eip);
  }

}

static Int4 CountNucleotides (CharPtr sequence)
{
  Int4    num = 0;
  CharPtr cp;
  
  if (sequence == NULL) return 0;	
  for (cp = sequence; *cp != 0; cp++)
  {
  	if (*cp != '-') 
  	{
  	  num++;
  	}
  }
  return num;
}

static void PrintAlignmentSummary (TAlignmentFilePtr afp, FILE *fp)
{
  Int4 index;

  if (fp == NULL) return;

  if (afp == NULL) {
    fprintf (fp, "Catastrophic failure during reading\n");
  } else {
    fprintf (fp, "Found %d sequences\n", afp->num_sequences);
    fprintf (fp, "Found %d organisms\n", afp->num_organisms);
    if (afp->num_sequences == afp->num_segments * afp->num_organisms)
    {
      for (index = 0; index < afp->num_sequences; index++)
      {
        fprintf (fp, "\t%s\t%d nucleotides\t", afp->ids [index],
                 CountNucleotides (afp->sequences[index]));
        if (index / afp->num_segments < afp->num_organisms) {
          fprintf (fp, "%s\n", afp->organisms [index / afp->num_segments]);
        } else {
          fprintf (fp, "No organism information\n");
        }
      }    	
    }
    else
    {
      for (index = 0; index < afp->num_sequences; index++)
      {
        fprintf (fp, "\t%s\t%d nucleotides\t", afp->ids [index], 
                 CountNucleotides (afp->sequences[index]));
        if (index < afp->num_organisms) {
          fprintf (fp, "%s\n", afp->organisms [index]);
        } else {
          fprintf (fp, "No organism information\n");
        }
      }
      while (index < afp->num_organisms) {
        fprintf (fp, "Unclaimed organism: %s\n", afp->organisms [index]);
        index++;
      }	
    }
  }
}


static void ReportPotentialDupIDs (TAlignmentFilePtr afp, FILE *fp)
{
  int     seq_index, k;
  int     curr_seg;
  int     num_sequences;
  Int4Ptr seq_len;
  BoolPtr may_be_dup;
  Int4    a, b;

  if (afp == NULL || afp->sequences == NULL || afp->num_sequences == 0) {
    return;
  }

  num_sequences = afp->num_sequences / afp->num_segments;

  may_be_dup = (BoolPtr) MemNew (sizeof (Boolean) * num_sequences);
  for (seq_index = 0; seq_index < num_sequences; seq_index++) {
    may_be_dup[seq_index] = FALSE;
  }

  seq_len = (Int4Ptr) MemNew (sizeof (Int4) * afp->num_sequences);
  for (seq_index = 0; seq_index < afp->num_sequences; seq_index++) {
    seq_len[seq_index] = StringLen (afp->sequences[seq_index]);
  }

  for (curr_seg = 0; curr_seg < afp->num_segments; curr_seg++) {
    for (seq_index = 0; seq_index < num_sequences - 1; seq_index++) {
      for (k = seq_index + 1; k < num_sequences; k++) {
        a = curr_seg * num_sequences + seq_index;
        b = curr_seg * num_sequences + k;
        if (seq_len[a] != seq_len[b]) {
          if (seq_len[a] % seq_len [b] == 0) {
            may_be_dup[seq_index] = TRUE;
          } else if (seq_len[b] % seq_len[a] == 0) {
            may_be_dup[k] = TRUE;
          }
        }
      }
    }
  }
  seq_len = MemFree (seq_len);

  for (seq_index = 0; seq_index < num_sequences; seq_index ++)
  {
    if (may_be_dup[seq_index]) {
      fprintf (fp, "Please check your file - %s may have been used as an ID for multiple sequences.\n", afp->ids[seq_index]);
    }
  }

  may_be_dup = MemFree (may_be_dup);
}


static Boolean DoSequenceLengthsMatch (TAlignmentFilePtr afp)
{
  int     seq_index;
  int     curr_seg;
  Int4Ptr seq_len;
  Boolean rval;

  if (afp == NULL || afp->sequences == NULL || afp->num_sequences == 0) {
    return TRUE;
  }

  seq_len = (Int4Ptr) MemNew (sizeof (Int4) * afp->num_segments);
  if (seq_len == NULL) return FALSE;
  for (seq_index = 0; seq_index < afp->num_segments; seq_index ++)
  {
    seq_len [seq_index] = StringLen (afp->sequences[seq_index]);
  }

  curr_seg = 0;
  rval = TRUE;
  for (seq_index = afp->num_segments; seq_index < afp->num_sequences && rval; seq_index++) {
    if (StringLen (afp->sequences[seq_index]) != seq_len[curr_seg]) {
      rval = FALSE;
    }
	curr_seg ++;
	if (curr_seg >= afp->num_segments) curr_seg = 0;
  }
  MemFree (seq_len);
  return rval;
}


extern void 
ProduceAlignmentNotes 
(TAlignmentFilePtr afp,
 TErrorInfoPtr error_list)
{
  Char         path [PATH_MAX];
  FILE         *fp;
  Boolean      ok_to_import = FALSE;

  TmpNam (path);
  fp = FileOpen (path, "wb");
  if (fp == NULL) return;


  if (afp != NULL && DoSequenceLengthsMatch (afp)) {
    ok_to_import = TRUE;
  }

  
  if (ok_to_import && error_list != NULL) {
    fprintf (fp, "Congratulations, you have successfully created a sequin file;\nhowever, I had trouble reading part of your file.\nPlease check your data carefully before submitting to be sure that all of your sequences\nwere included correctly.\nIf your file is incomplete, or contains incorrect sequences, please use the error report below\nto find the problem.\n");
  }
  ReportPotentialDupIDs (afp, fp);

  WalkErrorList (error_list, fp);
  PrintAlignmentSummary (afp, fp);

  FileClose (fp);
  LaunchGeneralTextViewer (path, "Alignment reading summary");
  FileRemove (path);
}

typedef struct alphabetformdata {
  FEATURE_FORM_BLOCK

  DialoG aln_settings;
  Handle obj;
  Char  path [PATH_MAX];
  FILE  *fp;
} AlphabetFormData, PNTR AlphabetFormPtr;

extern SeqEntryPtr 
SeqEntryFromAlignmentFile 
(FILE *fp,
 TSequenceInfoPtr sequence_info,
 Uint1   moltype,
 CharPtr no_org_err_msg)
{
  TErrorInfoPtr     error_list;
  ReadBufferData    rbd;
  TAlignmentFilePtr afp;
  SeqEntryPtr       sep = NULL;

  if (fp == NULL || sequence_info == NULL) return NULL;

  error_list = NULL;
  rbd.fp = fp;
  rbd.current_data = NULL;
  afp = ReadAlignmentFile ( AbstractReadFunction,
                            (Pointer) &rbd,
                            AbstractReportError,
                            (Pointer) &error_list,
                            sequence_info);

  ProduceAlignmentNotes (afp, error_list);
  ErrorInfoFree (error_list);
  if (afp != NULL) {
    if (afp->num_organisms == 0 && no_org_err_msg != NULL) {
      Message (MSG_ERROR, no_org_err_msg);
    } else if (afp->num_organisms != 0 && afp->num_organisms != afp->num_sequences && afp->num_organisms * afp->num_segments != afp->num_sequences) {
      Message (MSG_ERROR, "Number of organisms must match number of sequences!");
    } else if (! DoSequenceLengthsMatch (afp)) {
      Message (MSG_ERROR, "Your alignment is incorrectly formatted.  Sequence plus gaps should be the same length for all sequences.");
    } else {
      sep = MakeSequinDataFromAlignment (afp, moltype);
    }
  }

  AlignmentFileFree (afp);
  return sep;
}

static void DoReadAlignment (ButtoN b)
{
  AlphabetFormPtr  abc;
  TAlignmentFilePtr afp;
  TErrorInfoPtr     error_list;
  TSequenceInfoPtr  sequence_info;
  SeqEntryPtr      sep;
  Uint1            moltype;
  ReadBufferData   rbd;

  abc = GetObjectExtra (b);
  if (abc == NULL) return;
  Hide (abc->form);
  Update ();

  sequence_info = (TSequenceInfoPtr) DialogToPointer (abc->aln_settings);
  if (sequence_info == NULL) return;

  WatchCursor ();
  Update ();

  if (StringCmp (sequence_info->alphabet, protein_alphabet) == 0) {
    moltype = Seq_mol_aa;
  } else {
    moltype = Seq_mol_na;
  }

  error_list = NULL;
  rbd.fp = abc->fp;
  rbd.current_data = NULL;
  afp = ReadAlignmentFile ( AbstractReadFunction,
                            (Pointer) &rbd,
                            AbstractReportError,
                            (Pointer) &error_list,
                            sequence_info);
  if (afp != NULL) {
      sep = MakeSequinDataFromAlignment (afp, moltype);
      FinishAlignmentRead (abc->obj, sep, abc->path);
  }
  ProduceAlignmentNotes (afp, error_list);
  ErrorInfoFree (error_list);
  SequenceInfoFree (sequence_info);

  AlignmentFileFree (afp);
  FileClose (abc->fp);
  abc->fp = NULL;
  ArrowCursor ();
  Update ();
  Remove (abc->form);
}


/* Need cleanup for Alphabet Dialog to close File Pointer */
static void CleanupAlphabetDialog (GraphiC g, VoidPtr data)
{
    AlphabetFormPtr afp;

    afp = (AlphabetFormPtr) data;
    if (afp != NULL && afp->fp != NULL) {
        FileClose (afp->fp);
        afp->fp = NULL;
    }
    StdCleanupFormProc (g, data);
}

static void BuildGetAlphabetDialog (IteM i)
{
  BaseFormPtr        bfp;
  ButtoN             b;
  GrouP              c, h;
  AlphabetFormPtr    afp;
  WindoW             w;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;

  afp = (AlphabetFormPtr) MemNew (sizeof (AlphabetFormData));
  if (afp == NULL) return;
  afp->obj = i;

  if (! GetInputFileName (afp->path, sizeof (afp->path), NULL, "TEXT")) return;

  afp->fp = FileOpen (afp->path, "r");
  if (afp->fp == NULL)
  {
    Message (MSG_ERROR, "Unable to open file");
    return;
  }

  w = FixedWindow (-50, -33, -10, -10, "Alignment Alphabet", StdCloseWindowProc);
  SetObjectExtra (w, afp, CleanupAlphabetDialog);
  afp->form = (ForM) w;
  afp->formmessage = NULL;

  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  afp->aln_settings = AlnSettingsDlg (h, TRUE);  

  c = HiddenGroup (h, 4, 0, NULL);
  b = DefaultButton (c, "Accept", DoReadAlignment);
  SetObjectExtra (b, afp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) afp->aln_settings, (HANDLE) c, NULL);
  RealizeWindow (w);
  Show (w);
  Update ();
}

extern void ReadAlignment (IteM i)
{
  BuildGetAlphabetDialog (i);

}


typedef struct alnsettingsdlg {
  DIALOG_MESSAGE_BLOCK
  TexT missing;
  TexT match;
  TexT beginning_gap;
  TexT middle_gap;
  TexT end_gap;

  PopuP sequence_type;
} AlnSettingsDlgData, PNTR AlnSettingsDlgPtr;

static CharPtr aln_settings_help = "\
Beginning Gap: When some of the sequences in an \
alignment are shorter or longer than others, beginning \
gap characters are added to the beginning of the sequence \
to maintain the correct spacing.  These will not appear \
in your sequence file.\n\
Middle Gap: These characters are used to maintain the spacing \
inside an alignment.  These are not nucleotides and will \
not appear as part of your sequence file.\n\
End Gap: When some of the sequences in an alignment are shorter \
or longer than others, end gap characters are added to the end \
of the sequence to maintain the correct spacing.  These will \
not appear in your sequence file.\n\
Ambiguous/Unknown: These characters are used to represent \
indeterminate/ambiguous nucleotides.  These will appear in your \
sequence file as 'n'.\n\
Match: These characters are used to indicate positions where \
sequences are identical to the first sequence.  These will be \
replaced by the actual characters from the first sequence.";

const char *nucleotide_alphabet = "ABCDGHKMRSTUVWYabcdghkmrstuvwy";
const char *protein_alphabet = "ABCDEFGHIKLMPQRSTUVWXYZabcdefghiklmpqrstuvwxyz";

extern TSequenceInfoPtr GetDefaultSequenceInfo (void)
{
  TSequenceInfoPtr sequence_info = SequenceInfoNew();

  sequence_info->missing = MemFree (sequence_info->missing);
  sequence_info->missing = StringSave ("?Nn");

  sequence_info->beginning_gap = MemFree (sequence_info->beginning_gap);
  sequence_info->beginning_gap = StringSave ("-.Nn?");

  sequence_info->middle_gap = MemFree (sequence_info->middle_gap);
  sequence_info->middle_gap = StringSave ("-.");

  sequence_info->end_gap = MemFree (sequence_info->end_gap);
  sequence_info->end_gap = StringSave ("-.Nn?");

  sequence_info->match = MemFree (sequence_info->match);
  sequence_info->match = StringSave (":");

  sequence_info->alphabet = nucleotide_alphabet;

  return sequence_info;
}

static Pointer AlnSettingsDlgToData (DialoG d)
{
  AlnSettingsDlgPtr dlg;
  TSequenceInfoPtr  sequence_info;

  dlg = (AlnSettingsDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return NULL;

  sequence_info = SequenceInfoNew ();
  if (sequence_info == NULL) return NULL;

  sequence_info->missing = MemFree (sequence_info->missing);
  sequence_info->missing = SaveStringFromText (dlg->missing);

  sequence_info->beginning_gap = MemFree (sequence_info->beginning_gap);
  sequence_info->beginning_gap = SaveStringFromText (dlg->beginning_gap);

  sequence_info->middle_gap = MemFree (sequence_info->middle_gap);
  sequence_info->middle_gap = SaveStringFromText (dlg->middle_gap);

  sequence_info->end_gap = MemFree (sequence_info->end_gap);
  sequence_info->end_gap = SaveStringFromText (dlg->end_gap);

  sequence_info->match = MemFree (sequence_info->match);
  sequence_info->match = SaveStringFromText (dlg->match);

  if (dlg->sequence_type != NULL) 
  {
    if (GetValue (dlg->sequence_type) == 1) {
      sequence_info->alphabet = nucleotide_alphabet;
    } else {
      sequence_info->alphabet = protein_alphabet;
    }
  }
  else
  {
    sequence_info->alphabet = nucleotide_alphabet;
  }

  return sequence_info;
}


static void DataToAlnSettingsDlg (DialoG d, Pointer data)
{
  AlnSettingsDlgPtr dlg;
  TSequenceInfoPtr  sequence_info;

  dlg = (AlnSettingsDlgPtr) GetObjectExtra (d);
  if (dlg == NULL) return;

  sequence_info = (TSequenceInfoPtr) data;

  if (sequence_info == NULL)
  {
    SetTitle (dlg->missing, "?Nn");
    SetTitle (dlg->beginning_gap, "-.Nn?");
    SetTitle (dlg->middle_gap, "-.");
    SetTitle (dlg->end_gap, "-.Nn?");
    SetTitle (dlg->match, ":");
    if (dlg->sequence_type != NULL)
    {
      SetValue (dlg->sequence_type, 1);
    }
  }
  else
  {
    SetTitle (dlg->missing, sequence_info->missing);
    SetTitle (dlg->beginning_gap, sequence_info->beginning_gap);
    SetTitle (dlg->middle_gap, sequence_info->middle_gap);
    SetTitle (dlg->end_gap, sequence_info->end_gap);
    SetTitle (dlg->match, sequence_info->match);

    if (dlg->sequence_type != NULL) 
    {
      if (StringCmp (sequence_info->alphabet, protein_alphabet) == 0) 
      {
        SetValue (dlg->sequence_type, 2);
      }
      else
      {
        SetValue (dlg->sequence_type, 1);
      }
    }
  }
}


static ValNodePtr TestAlnSettingsDlg (DialoG d)
{
  ValNodePtr        err_list = NULL;
  TSequenceInfoPtr  sequence_info;
  CharPtr           cp;
  CharPtr           fmt = "Character %c cannot appear in both %s and %s.";
  CharPtr           err_str;
  CharPtr           missing_name = "Ambiguous/Unknown";
  CharPtr           middle_gap_name = "Middle Gap";
  CharPtr           match_name = "Match";

  sequence_info = DialogToPointer (d);
  if (sequence_info == NULL) return NULL;

  /* missing and match cannot appear in middle gap list */
  cp = sequence_info->missing;
  while (cp != NULL && *cp != 0)
  {
    if (StringChr (sequence_info->middle_gap, *cp)) 
    {
      err_str = (CharPtr) MemNew (sizeof (Char) * (StringLen (fmt) 
                                                   + StringLen (missing_name)
                                                   + StringLen (middle_gap_name)));
      sprintf (err_str, fmt, *cp, missing_name, middle_gap_name);
      ValNodeAddPointer (&err_list, 0, err_str);
    }
    cp++;
  }

  cp = sequence_info->match;
  while (cp != NULL && *cp != 0)
  {
    if (StringChr (sequence_info->middle_gap, *cp)) 
    {
      err_str = (CharPtr) MemNew (sizeof (Char) * (StringLen (fmt) 
                                                   + StringLen (match_name)
                                                   + StringLen (middle_gap_name)));
      sprintf (err_str, fmt, *cp, match_name, middle_gap_name);
      ValNodeAddPointer (&err_list, 0, err_str);
    }
    cp++;
  }

  /* missing and match cannot share characters */
  cp = sequence_info->missing;
  while (cp != NULL && *cp != 0)
  {
    if (StringChr (sequence_info->match, *cp)) 
    {
      err_str = (CharPtr) MemNew (sizeof (Char) * (StringLen (fmt) 
                                                   + StringLen (missing_name)
                                                   + StringLen (match_name)));
      sprintf (err_str, fmt, *cp, missing_name, match_name);
      ValNodeAddPointer (&err_list, 0, err_str);
    }
    cp++;
  }

  return err_list;
}


extern DialoG AlnSettingsDlg (GrouP h, Boolean allow_sequence_type)
{
  AlnSettingsDlgPtr dlg;
  GrouP             p, g, p_msg;
  
  dlg = (AlnSettingsDlgPtr) MemNew (sizeof (AlnSettingsDlgData));
  p = HiddenGroup (h, -1, 0, NULL);
  SetObjectExtra (p, dlg, StdCleanupExtraProc);
  SetGroupSpacing (p, 10, 10);

  dlg->dialog = (DialoG) p;
  dlg->todialog = DataToAlnSettingsDlg;
  dlg->fromdialog = AlnSettingsDlgToData;
  dlg->testdialog = TestAlnSettingsDlg;

  g = HiddenGroup (p, 2, 0, NULL);
  StaticPrompt (g, "Ambiguous/Unknown", 0, dialogTextHeight, programFont, 'c');
  dlg->missing = DialogText (g, "?Nn", 5, NULL);
  StaticPrompt (g, "Match", 0, dialogTextHeight, programFont, 'c');
  dlg->match = DialogText (g, ".", 5, NULL);
  StaticPrompt (g, "Beginning Gap", 0, dialogTextHeight, programFont, 'c');
  dlg->beginning_gap = DialogText (g, "-.?nN", 5, NULL);
  StaticPrompt (g, "Middle Gap", 0, dialogTextHeight, programFont, 'c');
  dlg->middle_gap = DialogText (g, "-", 5, NULL);
  StaticPrompt (g, "End Gap", 0, dialogTextHeight, programFont, 'c');
  dlg->end_gap = DialogText (g, "-.?nN", 5, NULL);
  if (allow_sequence_type) {
    StaticPrompt (g, "Sequence Type", 0, dialogTextHeight, programFont, 'c');
    dlg->sequence_type = PopupList (g, TRUE, NULL);
    PopupItem (dlg->sequence_type, "Nucleotide");
    PopupItem (dlg->sequence_type, "Protein");
    SetValue (dlg->sequence_type, 1);
  }
  
  p_msg = MultiLinePrompt (p, aln_settings_help, 750, systemFont);

  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) p_msg, NULL);
  
  return (DialoG) p;
}


static void ReadNewAsnProc (IteM i)

{
  CommonReadNewAsnProc ((Handle) i, FALSE, FALSE);
}

extern Boolean LIBCALLBACK SequinOpenMimeFile (CharPtr filename)

{
  Boolean  rsult = FALSE;

#ifdef WIN_MAC
  SafeHide (startupForm); /* if just Hide, triggers MacDeactProc, does not reactivate */
  Update ();
#endif
  rsult = DoReadAnythingLoop (NULL, filename, filename, FALSE, FALSE, FALSE, TRUE, FALSE);
  if (! rsult) {
/*
#ifndef WIN16
    if (BiostrucAvail () && OpenMimeFileWithDeletion (filename, FALSE)) {
      Handle   www_cn3d;
      www_cn3d = Cn3DWin_Entrez(NULL, useEntrez);
      Show (www_cn3d);
      Select (www_cn3d);
    } else {
      Show (startupForm);
      Select (startupForm);
    }
#else
    Show (startupForm);
    Select (startupForm);
#endif
*/
    Show (startupForm);
    Select (startupForm);
  }
#ifdef OS_MAC
  /* the Web browsers on other platforms get upset if you delete the file, */
  /* but apparently on the Mac you must delete it - comment from cn3dmain.c */
  /* but here it's for response to apple events, so should not free this file */
  /* FileRemove (filename); */
#endif
  return rsult;
}

extern Boolean LIBCALLBACK SequinOpenResultFile (CharPtr filename)

{
  return DoReadAnythingLoop (NULL, filename, NULL, FALSE, FALSE, FALSE, TRUE, FALSE);
}

extern Boolean LIBCALLBACK SequinHandleNetResults (CharPtr filename)

{
  return DoReadAnythingLoop (NULL, filename, NULL, FALSE, FALSE, TRUE, TRUE, TRUE);
}

#ifdef WIN_MAC
static void MacReadNewAsnProc (IteM i)

{
  if (initialFormsActive) {
    CommonReadNewAsnProc ((Handle) i, TRUE, FALSE);
  } else {
    CommonReadNewAsnProc ((Handle) i, FALSE, FALSE);
  }
}
#endif

static Boolean CountAlignmentsCallback (GatherContextPtr gcp)

{
  Int4Ptr  rsultptr;

  if (gcp == NULL) return TRUE;

  rsultptr = (Int4Ptr) gcp->userdata;
  if (rsultptr == NULL ) return TRUE;

  switch (gcp->thistype) {
    case OBJ_SEQALIGN :
    case OBJ_SEQHIST_ALIGN :
      (*rsultptr)++;
      return TRUE;
    default :
      break;
  }
  return TRUE;
}

static Int4 LIBCALL CountSeqEntryAligns (Uint2 entityID, SeqEntryPtr sep)

{
  GatherScope  gs;
  Int4         rsult;

  rsult = 0;
  if (entityID == 0 || sep == NULL) return 0;
  MemSet ((Pointer) (&gs), 0, sizeof (GatherScope));
  gs.seglevels = 1;
  MemSet((Pointer) (gs.ignore), (int) (TRUE), (size_t) (OBJ_MAX * sizeof (Boolean)));
  gs.ignore[OBJ_BIOSEQ] = FALSE;
  gs.ignore[OBJ_BIOSEQ_SEG] = FALSE;
  gs.ignore[OBJ_SEQALIGN] = FALSE;
  gs.ignore[OBJ_SEQANNOT] = FALSE;
  gs.ignore[OBJ_SEQHIST] = FALSE;
  gs.ignore[OBJ_SEQHIST_ALIGN] = FALSE;
  gs.scope = sep;
  GatherEntity (entityID, (Pointer) (&rsult), CountAlignmentsCallback, &gs);
  return rsult;
}

static BioseqSetPtr GetAlignmentSetForBioseq (BioseqPtr bsp, Uint2 entityID)
{
  BioseqSetPtr bssp = NULL;
  SeqEntryPtr  sep, sep_top;
  Uint2        parenttype;
  Pointer      parentptr;
  
  if (bsp == NULL)
  {
    return NULL;
  }
  
  sep_top = GetTopSeqEntryForEntityID (entityID);
  
  sep = GetBestTopParentForData (entityID, bsp); 
  GetSeqEntryParent (sep, &parentptr, &parenttype);
  if (parenttype == OBJ_BIOSEQSET)
  {
    bssp = (BioseqSetPtr) parentptr;
  }
  else if (sep != NULL && IS_Bioseq_set (sep))
  {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
  }
  else if (sep_top != NULL && IS_Bioseq_set (sep_top))
  {
    bssp = (BioseqSetPtr) sep_top->data.ptrvalue;
  }
  return bssp;
}

static void EnableEditAlignItem (BaseFormPtr bfp)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  IteM          editalign;
  IteM          editdup;
  IteM          editfeatprop;
  Int2          mssgalign;
  Int2          mssgdup;
  Int2          mssgfeatprop;
  Int4          num;
  SeqEntryPtr   sep;
  SelStructPtr  sel;

  if (bfp == NULL && bfp->input_entityID != 0) return;
  mssgalign = RegisterFormMenuItemName ("SequinEditAlignmentItem");
  mssgdup = RegisterFormMenuItemName ("SequinDuplicateItem");
  mssgfeatprop = RegisterFormMenuItemName ("SequinFeaturePropagate");
  editalign = FindFormMenuItem (bfp, mssgalign);
  editdup = FindFormMenuItem (bfp, mssgdup);
  editfeatprop = FindFormMenuItem (bfp, mssgfeatprop);
  sel = ObjMgrGetSelected ();
  sep = NULL;
  /*
  bsp =  GetBioseqGivenIDs (bfp->input_entityID, bfp->input_itemID, bfp->input_itemtype);
  if (bsp != NULL) {
    sep = GetBestTopParentForData (bfp->input_entityID, bsp);
  } else {
    sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  }
  */
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  num = CountSeqEntryAligns (bfp->input_entityID, sep);
  if (sep != NULL && num > 0) {
    if (sel != NULL && sel->entityID == bfp->input_entityID &&
        (sel->itemtype == OBJ_SEQALIGN || sel->itemtype == OBJ_SEQHIST_ALIGN)) {
      Enable (editalign);
    } else /* if (sel == NULL) */ {
      if (FindSeqAlignInSeqEntry (sep, OBJ_SEQALIGN) != NULL) {
        Enable (editalign);
      } else {
        Disable (editalign);
      }
    } /* else {
      Disable (editalign);
      Disable (editupwthaln);
    } */
  } else {
    Disable (editalign);
  }
  if (sel != NULL &&
      (sel->itemtype == OBJ_SEQDESC || sel->itemtype == OBJ_SEQFEAT)) {
    Enable (editdup);
  } else {
    Disable (editdup);
  }
  bsp = GetBioseqGivenIDs (bfp->input_entityID, bfp->input_itemID, bfp->input_itemtype);
  if (num > 0 && bsp != NULL && sep != NULL && IS_Bioseq_set (sep)) {
    bssp = GetAlignmentSetForBioseq (bsp, bfp->input_entityID);
    if (bssp != NULL) {
      if (IsPopPhyEtcSet (bssp->_class)) {
        Enable (editfeatprop);
      } else if (bssp->_class == 7 && indexerVersion) {
        Enable (editfeatprop);
      } else {
        Disable (editfeatprop);
      }
    } else {
      Disable (editfeatprop);
    }
  } else {
    Disable (editfeatprop);
  }
}

static void EnableEditSeqAlignAndSubItems (BaseFormPtr bfp)

{
  BioseqPtr      bsp;
  BioseqSetPtr   bssp;
  IteM           editseq;
  IteM           editsub;
  MenU           editupd;
  MenU           editupd_idx = NULL;
  MenU           editext;
  MenU           editadd;
  Int2           mssgadd;
  Int2           mssgseq;
  Int2           mssgsub;
  Int2           mssgupd;
  Int2           mssgupd_idx = 0;
  Int2           mssgext;
  ObjMgrDataPtr  omdp;
  SeqEntryPtr    sep;

  if (bfp == NULL && bfp->input_entityID != 0) return;
  mssgseq = RegisterFormMenuItemName ("SequinEditSequenceItem");
  mssgsub = RegisterFormMenuItemName ("SequinEditSubmitterItem");
  mssgupd = RegisterFormMenuItemName ("SequinUpdateSeqSubmenu");
  if (indexerVersion) {
    mssgupd_idx = RegisterFormMenuItemName ("SequinUpdateSeqSubmenuIndexer");
  }
  mssgext = RegisterFormMenuItemName ("SequinExtendSeqSubmenu");
  mssgadd = RegisterFormMenuItemName ("SequinAddSeqSubmenu");
  editseq = FindFormMenuItem (bfp, mssgseq);
  editsub = FindFormMenuItem (bfp, mssgsub);
  editupd = (MenU) FindFormMenuItem (bfp, mssgupd);
  if (indexerVersion) {
    editupd_idx = (MenU) FindFormMenuItem (bfp, mssgupd_idx);
  }
  editext = (MenU) FindFormMenuItem (bfp, mssgext);
  editadd = (MenU) FindFormMenuItem (bfp, mssgadd);
  bsp = GetBioseqGivenIDs (bfp->input_entityID, bfp->input_itemID, bfp->input_itemtype);
  if (bsp != NULL) {
    Enable (editseq);
    Enable (editupd);
    if (indexerVersion) {
      Enable (editupd_idx);
    }
    Enable (editext);
    sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
    if (sep != NULL && IS_Bioseq_set (sep)) {
      bssp = (BioseqSetPtr) sep->data.ptrvalue;
      if (bssp != NULL && (bssp->_class == 7 || (IsPopPhyEtcSet (bssp->_class)))) {
        Enable (editadd);
      } else {
        Disable (editadd);
      }
    } else {
      Disable (editadd);
    }
#ifdef WIN_MAC
    Enable (featPropItem);
#endif
  } else {
    Disable (editseq);
    Disable (editupd);
    if (indexerVersion) {
      Disable (editupd_idx);
    }
    Disable (editext);
    Disable (editadd);
#ifdef WIN_MAC
    Disable (featPropItem);
#endif
  }
  omdp = ObjMgrGetData (bfp->input_entityID);
  if (omdp != NULL && omdp->datatype == OBJ_SEQSUB) {
    Enable (editsub);
#ifdef WIN_MAC
    Enable (prepareItem);
    Enable (submitItem);
#endif
  } else {
    Disable (editsub);
#ifdef WIN_MAC
    Disable (prepareItem);
    Disable (submitItem);
#endif
  }
  EnableEditAlignItem (bfp);
}

#ifdef WIN_MAC
static void MedlineViewFormActivated (WindoW w)

{
  currentFormDataPtr = (VoidPtr) GetObjectExtra (w);
  initialFormsActive = FALSE;
  RepeatProcOnHandles (Enable,
                   (HANDLE) closeItem,
                   (HANDLE) duplicateViewItem,
                   (HANDLE) exportItem,
                   (HANDLE) printItem,
                   (HANDLE) copyItem,
                   (HANDLE) displayfontItem,
                   NULL);
}

static void BioseqViewFormActivated (WindoW w)

{
  bioseqViewUp = TRUE;
  currentFormDataPtr = (VoidPtr) GetObjectExtra (w);
  initialFormsActive = FALSE;
  /* SafeSetTitle (importItem, "Import Nucleotide FASTA..."); */
  RepeatProcOnHandles (Enable,
                   (HANDLE) openItem,
                   (HANDLE) closeItem,
                   (HANDLE) saveItem,
                   (HANDLE) saveAsItem,
                   (HANDLE) copyItem,
                   (HANDLE) deleteItem,
                   (HANDLE) findItem,
                   (HANDLE) findFFItem,
                   (HANDLE) findGeneItem,
                   (HANDLE) findProtItem,
                   (HANDLE) findPosItem,
                   (HANDLE) exportItem,
                   (HANDLE) duplicateViewItem,
                   (HANDLE) restoreItem,
                   (HANDLE) printItem,
                   (HANDLE) orfItem,
                   (HANDLE) targetItem,
                   (HANDLE) newDescMenu,
                   (HANDLE) newFeatMenu,
                   (HANDLE) advTableMenu,
                   (HANDLE) sucItem,
                   (HANDLE) newPubMenu,
                   (HANDLE) batchApplyMenu,
                   (HANDLE) batchEditMenu,
                   (HANDLE) prepareItem,
                   (HANDLE) edithistoryitem,
                   NULL);
  Enable (validateItem);
  Enable (validateMenu);
  Enable (vecscreenMenu);
  Enable (aluItem);
  Enable (submitItem);
  Enable (specialMenu);
  Enable (projectsMenu);
  Enable (analysisMenu);
  Enable (vectorScreenItem);
  Enable (cddBlastItem);
  Enable (cddSearchMenu);
  Enable (cddSearchItem);
  Enable (parseFileItem);
  Enable (spellItem);
  Enable (addSecondaryItem);
  Enable (displayfontItem);
  Disable (importItem);
  EnableFeaturesPerTarget ((BaseFormPtr) currentFormDataPtr);
  EnableAnalysisItems ((BaseFormPtr) currentFormDataPtr, FALSE);
  EnableEditSeqAlignAndSubItems ((BaseFormPtr) currentFormDataPtr);
}

/*
static void TermSelectionActivateProc (WindoW w)

{
  termListUp = TRUE;
  currentFormDataPtr = (VoidPtr) GetObjectExtra (w);
  initialFormsActive = FALSE;
  if (UsingTextQuery ((ForM) w)) {
    SafeSetValue (queryChoice, 2);
  } else {
    SafeSetValue (queryChoice, 1);
  }
  RepeatProcOnHandles (Enable,
                   (HANDLE) importItem,
                   (HANDLE) exportItem,
                   (HANDLE) cutItem,
                   (HANDLE) copyItem,
                   (HANDLE) pasteItem,
                   (HANDLE) deleteItem,
                   (HANDLE) preferencesItem,
                   (HANDLE) queryChoice,
                   (HANDLE) clearUnusedItem,
                   NULL);
  Enable (loadUidItem);
  Enable (saveUidItem);
}

static void DocumentSummaryActivateProc (WindoW w)

{
  docSumUp = TRUE;
  currentFormDataPtr = (VoidPtr) GetObjectExtra (w);
  initialFormsActive = FALSE;
  if (UsingDelayedNeighbor ((ForM) w)) {
    SafeSetValue (neighborChoice, 2);
  } else {
    SafeSetValue (neighborChoice, 1);
  }
  RepeatProcOnHandles (Enable,
                   (HANDLE) saveAsItem,
                   (HANDLE) closeItem,
                   (HANDLE) importItem,
                   (HANDLE) exportItem,
                   (HANDLE) printItem,
                   (HANDLE) copyItem,
                   (HANDLE) preferencesItem,
                   (HANDLE) neighborChoice,
                   (HANDLE) docsumfontItem,
                   (HANDLE) analysisMenu,
                   NULL);
  Enable (loadUidItem);
  Enable (saveUidItem);
  EnableAnalysisItems ((BaseFormPtr) currentFormDataPtr, TRUE);
}
*/

static void SeqEditFormActivated (WindoW w)

{
  currentFormDataPtr = (VoidPtr) GetObjectExtra (w);
  initialFormsActive = FALSE;
  RepeatProcOnHandles (Enable,
                   (HANDLE) openItem,
                   (HANDLE) closeItem,
                   (HANDLE) cutItem,
                   (HANDLE) copyItem,
                   (HANDLE) pasteItem,
                   (HANDLE) deleteItem,
                   (HANDLE) importItem,
                   (HANDLE) exportItem,
                   NULL);
}

static void StdEditorFormActivated (WindoW w)

{
  currentFormDataPtr = (VoidPtr) GetObjectExtra (w);
  initialFormsActive = FALSE;
  RepeatProcOnHandles (Enable,
                   (HANDLE) openItem,
                   (HANDLE) closeItem,
                   (HANDLE) cutItem,
                   (HANDLE) copyItem,
                   (HANDLE) pasteItem,
                   (HANDLE) deleteItem,
                   (HANDLE) importItem,
                   (HANDLE) exportItem,
                   NULL);
}

static void StdValidatorFormActivated (WindoW w)

{
  currentFormDataPtr = (VoidPtr) GetObjectExtra (w);
  initialFormsActive = FALSE;
  RepeatProcOnHandles (Enable,
                   (HANDLE) openItem,
                   (HANDLE) closeItem,
                   (HANDLE) cutItem,
                   (HANDLE) copyItem,
                   (HANDLE) pasteItem,
                   (HANDLE) deleteItem,
                   (HANDLE) importItem,
                   (HANDLE) exportItem,
                   (HANDLE) printItem,
                   NULL);
}

static void TextViewProcFormActivated (WindoW w)

{
  currentFormDataPtr = (VoidPtr) GetObjectExtra (w);
  initialFormsActive = FALSE;
  RepeatProcOnHandles (Enable,
                   (HANDLE) openItem,
                   (HANDLE) closeItem,
                   (HANDLE) cutItem,
                   (HANDLE) copyItem,
                   (HANDLE) pasteItem,
                   (HANDLE) deleteItem,
                   (HANDLE) importItem,
                   (HANDLE) exportItem,
                   (HANDLE) printItem,
                   NULL);
}

static void ConfigFormActivated (WindoW w)

{
  currentFormDataPtr = (VoidPtr) GetObjectExtra (w);
  initialFormsActive = FALSE;
  RepeatProcOnHandles (Enable,
                   (HANDLE) cutItem,
                   (HANDLE) copyItem,
                   (HANDLE) pasteItem,
                   (HANDLE) deleteItem,
                   NULL);
}

static void StartupActivateProc (WindoW w)

{
  currentFormDataPtr = (VoidPtr) GetObjectExtra (w);
  initialFormsActive = FALSE;
  Enable (openItem);
}

static void FormatActivateProc (WindoW w)

{
  currentFormDataPtr = (VoidPtr) GetObjectExtra (w);
  initialFormsActive = TRUE;
  Enable (openItem);
}

static void SubmitBlockActivateProc (WindoW w)

{
  currentFormDataPtr = (VoidPtr) GetObjectExtra (w);
  initialFormsActive = TRUE;
  RepeatProcOnHandles (Enable,
                   (HANDLE) openItem,
                   (HANDLE) cutItem,
                   (HANDLE) copyItem,
                   (HANDLE) pasteItem,
                   (HANDLE) deleteItem,
                   (HANDLE) importItem,
                   (HANDLE) exportItem,
                   NULL);
}

static void GenomeFormActivateProc (WindoW w)

{
  currentFormDataPtr = (VoidPtr) GetObjectExtra (w);
  initialFormsActive = TRUE;
  RepeatProcOnHandles (Enable,
                   (HANDLE) openItem,
                   (HANDLE) cutItem,
                   (HANDLE) copyItem,
                   (HANDLE) pasteItem,
                   (HANDLE) deleteItem,
                   (HANDLE) importItem,
                   (HANDLE) exportItem,
                   NULL);
}

static void OrgAndSeqsActivateProc (WindoW w)

{
  currentFormDataPtr = (VoidPtr) GetObjectExtra (w);
  initialFormsActive = TRUE;
  RepeatProcOnHandles (Enable,
                   (HANDLE) openItem,
                   (HANDLE) cutItem,
                   (HANDLE) copyItem,
                   (HANDLE) pasteItem,
                   (HANDLE) deleteItem,
                   (HANDLE) importItem,
                   (HANDLE) exportItem,
                   NULL);
}

static void HelpActivateProc (WindoW w)

{
  currentFormDataPtr = (VoidPtr) GetObjectExtra (w);
  RepeatProcOnHandles (Enable,
                   (HANDLE) openItem,
                   (HANDLE) closeItem,
                   (HANDLE) cutItem,
                   (HANDLE) copyItem,
                   (HANDLE) pasteItem,
                   (HANDLE) deleteItem,
                   (HANDLE) exportItem,
                   (HANDLE) printItem,
                   NULL);
}

static void MacDeactProc (WindoW w)

{
  termListUp = FALSE;
  docSumUp = FALSE;
  bioseqViewUp = FALSE;
  currentFormDataPtr = NULL;
  initialFormsActive = FALSE;
  Enable (openItem);
  SafeSetTitle (importItem, "Import...");
  SafeSetTitle (exportItem, "Export...");
  SafeSetValue (queryChoice, 0);
  SafeSetValue (neighborChoice, 0);
  RepeatProcOnHandles (Disable,
                   (HANDLE) saveItem,
                   (HANDLE) saveAsItem,
                   (HANDLE) closeItem,
                   (HANDLE) cutItem,
                   (HANDLE) copyItem,
                   (HANDLE) pasteItem,
                   (HANDLE) deleteItem,
                   (HANDLE) duplicateViewItem,
                   (HANDLE) importItem,
                   (HANDLE) exportItem,
                   (HANDLE) findItem,
                   (HANDLE) findFFItem,
                   (HANDLE) findGeneItem,
                   (HANDLE) findProtItem,
                   (HANDLE) findPosItem,
                   (HANDLE) orfItem,
                   (HANDLE) targetItem,
                   (HANDLE) newDescMenu,
                   (HANDLE) newFeatMenu,
                   (HANDLE) advTableMenu,
                   (HANDLE) sucItem,
                   (HANDLE) newPubMenu,
                   (HANDLE) batchApplyMenu,
                   (HANDLE) batchEditMenu,
                   (HANDLE) prepareItem,
                   (HANDLE) editsequenceitem,
                   (HANDLE) editseqalignitem,
                   (HANDLE) editseqsubitem,
                   (HANDLE) edithistoryitem,
                   (HANDLE) updateSeqMenu,
                   (HANDLE) extendSeqMenu,
                   (HANDLE) addSeqMenu,
                   (HANDLE) featPropItem,
                   (HANDLE) updalignitem,
                   (HANDLE) updalignidxitem,
                   NULL);
  Disable (validateItem);
  Disable (validateMenu);
  Disable (vecscreenMenu);
  Disable (aluItem);
  Disable (submitItem);
  Disable (vectorScreenItem);
  Disable (cddBlastItem);
  Disable (cddSearchMenu);
  Disable (cddSearchItem);
  Disable (parseFileItem);
  Disable (spellItem);
  Disable (docsumfontItem);
  Disable (queryChoice);
  Disable (clearUnusedItem);
  Disable (neighborChoice);
  Disable (legendItem);
  Disable (duplicateItem);
  Disable (addSecondaryItem);
  Disable (displayfontItem);
  Disable (printItem);
  Disable (restoreItem);
  Disable (specialMenu);
  Disable (projectsMenu);
  Disable (analysisMenu);
  Disable (loadUidItem);
  Disable (saveUidItem);
  if (indexerVersion) {
    Disable (updateSeqMenuIndexer);
  }
}
#endif

#ifndef WIN_MAC
#define MedlineViewFormActivated NULL
#define BioseqViewFormActivated NULL
#define TermSelectionActivateProc NULL
#define DocumentSummaryActivateProc NULL
#define SeqEditFormActivated NULL
#define StdEditorFormActivated NULL
#define StdValidatorFormActivated NULL
#define TextViewProcFormActivated NULL
#define ConfigFormActivated NULL
#define StartupActivateProc NULL
#define FormatActivateProc NULL
#define SubmitBlockActivateProc NULL
#define OrgAndSeqsActivateProc NULL
#define GenomeFormActivateProc NULL
#define HelpActivateProc NULL
#endif

static void HideHelpForm (ButtoN b)

{
  Hide (ParentWindow (b));
}

static void DisplayHelpFormProc (IteM i)

{
  if (helpForm == NULL) {
    WatchCursor ();
    helpForm = CreateHelpForm (-95, -5, "Sequin Help", "sequin.hlp",
                               HideHelpForm, HelpActivateProc);
    ArrowCursor ();
  }
  if (helpForm != NULL) {
    Show (helpForm);
    Select (helpForm);
  }
}

/*#ifdef USE_MEDARCH*/

static Boolean debug_fix_pub_equiv = FALSE;
static Boolean debug_fix_pub_set = FALSE;

static Boolean log_mla_asn = FALSE;
static Boolean log_mla_set = FALSE;

static ValNodePtr LookupAnArticleFuncNew (ValNodePtr oldpep, BoolPtr success)

{
  CitArtPtr      cap = NULL;
  ArticleIdPtr   ids;
  MlaBackPtr     mbp;
  MlaRequestPtr  mrp;
  Int4           pmid = 0;
  ValNodePtr     pub = NULL;
  ValNodePtr     vnp;
#ifdef OS_UNIX
  CharPtr        str;
#endif

  if (success != NULL) {
    *success = FALSE;
  }
  if (oldpep == NULL) return NULL;

#ifdef OS_UNIX
  if (! log_mla_set) {
    str = (CharPtr) getenv ("LOG_MLA_ASN");
    if (StringDoesHaveText (str)) {
      if (StringICmp (str, "TRUE") == 0) {
        log_mla_asn = TRUE;
      }
    }
    log_mla_set = TRUE;
  }
#endif

  for (vnp = oldpep; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == PUB_Article) {
      cap = (CitArtPtr) vnp->data.ptrvalue;
    } else if (vnp->choice == PUB_PMid) {
      pmid = (Int4) vnp->data.intvalue;
    }
  }

  if (cap != NULL) {
    mrp = Mla2CreateCitArtMatchRequest (cap);
    if (mrp != NULL) {
      if (log_mla_asn) {
        LaunchAsnTextViewer ((Pointer) mrp, (AsnWriteFunc) MlaRequestAsnWrite, "citart match request");
      }
      mbp = Mla2SynchronousQuery (mrp);
      mrp = MlaRequestFree (mrp);
      if (mbp != NULL) {
        if (log_mla_asn) {
          LaunchAsnTextViewer ((Pointer) mbp, (AsnWriteFunc) MlaBackAsnWrite, "citart match result");
        }
        pmid = Mla2ExtractCitMatchReply (mbp);
        mbp = MlaBackFree (mbp);
      }
    }
  }

  if (pmid > 0) {
    mrp = Mla2CreatePubFetchRequest (pmid);
    if (mrp != NULL) {
      if (log_mla_asn) {
        LaunchAsnTextViewer ((Pointer) mrp, (AsnWriteFunc) MlaRequestAsnWrite, "get pub request");
      }
      mbp = Mla2SynchronousQuery (mrp);
      mrp = MlaRequestFree (mrp);
      if (mbp != NULL) {
        if (log_mla_asn) {
          LaunchAsnTextViewer ((Pointer) mbp, (AsnWriteFunc) MlaBackAsnWrite, "get pub result");
        }
        cap = Mla2ExtractPubFetchReply (mbp);
        if (cap != NULL) {
          if (log_mla_asn) {
            LaunchAsnTextViewer ((Pointer) cap, (AsnWriteFunc) CitArtAsnWrite, "before");
          }
          ChangeCitArtMLAuthorsToSTD (cap);
          if (log_mla_asn) {
            LaunchAsnTextViewer ((Pointer) cap, (AsnWriteFunc) CitArtAsnWrite, "after");
          }
          for (ids = cap->ids; ids != NULL; ids = ids->next) {
            if (ids->choice != ARTICLEID_PUBMED) continue;
            if (ids->data.intvalue != pmid) {
              Message (MSG_POSTERR, "CitArt ID %ld does not match PMID %ld",
                       (long) ids->data.intvalue, (long) pmid);
            }
          }
          pub = ValNodeAddPointer (NULL, PUB_Article, (Pointer) cap);
          ValNodeAddInt (&pub, PUB_PMid, pmid);
          if (success != NULL) {
            *success = TRUE;
          }
        }
        mbp = MlaBackFree (mbp);
      }
    }
  }

  return pub;
}

static ValNodePtr LookupAnArticleFunc (ValNodePtr oldpep, BoolPtr success)

{
  FindPubOption  fpo;
  MonitorPtr     mon;
  Int4           muid = 0;
  ValNodePtr     pep;
  Int4           pmid = 0;
  ValNodePtr     pub;
  ValNodePtr     vnp;
#ifdef OS_UNIX
  AsnIoPtr       aip;
  CharPtr        str;
#endif

  pub = NULL;
  if (oldpep != NULL) {
    WatchCursor ();
    mon = MonitorStrNewEx ("Lookup Article", 40, FALSE);
    MonitorStrValue (mon, "Connecting to MedArch");
    Update ();
    if (MedArchInit ()) {
#ifdef OS_UNIX
      if (! debug_fix_pub_set) {
        str = (CharPtr) getenv ("DEBUG_FIX_PUB_EQUIV");
        if (StringDoesHaveText (str)) {
          if (StringICmp (str, "TRUE") == 0) {
            debug_fix_pub_equiv = TRUE;
          }
        }
        debug_fix_pub_set = TRUE;
      }
#endif
      pep = AsnIoMemCopy (oldpep, (AsnReadFunc) PubEquivAsnRead,
                          (AsnWriteFunc) PubEquivAsnWrite);
      fpo.always_look = TRUE;
      fpo.replace_cit = TRUE;
      fpo.lookups_attempted = 0;
      fpo.lookups_succeeded = 0;
      fpo.fetches_attempted = 0;
      fpo.fetches_succeeded = 0;
      MonitorStrValue (mon, "Performing Lookup");
      Update ();
#ifdef OS_UNIX
      if (debug_fix_pub_equiv) {
        if (pep != NULL) {
          aip = AsnIoOpen ("origpubequiv.txt", "w");
          if (aip != NULL) {
            PubEquivAsnWrite (pep, aip, NULL);
            AsnIoClose (aip);
          }
        }
      }
#endif
      pub = FixPubEquiv (pep, &fpo);
      if (fpo.fetches_succeeded && success != NULL) {
        *success = TRUE;
      }
#ifdef OS_UNIX
      if (debug_fix_pub_equiv) {
        if (pub != NULL) {
          aip = AsnIoOpen ("fixpubequiv.txt", "w");
          if (aip != NULL) {
            PubEquivAsnWrite (pub, aip, NULL);
            AsnIoClose (aip);
          }
        }
      }
#endif
      if (! fpo.fetches_succeeded) {
        ErrShow ();
        Update ();
      } else {
        for (vnp = pub; vnp != NULL; vnp = vnp->next) {
          if (vnp->choice == PUB_Muid) {
            muid = vnp->data.intvalue;
          } else if (vnp->choice == PUB_PMid) {
            pmid = vnp->data.intvalue;
          }
        }
        if (muid != 0 && pmid == 0) {
          pmid = MedArchMu2Pm (muid);
          if (pmid != 0) {
            ValNodeAddInt (&pub, PUB_PMid, pmid);
          }
        } else if (pmid != 0 && muid == 0) {
          muid = MedArchPm2Mu (pmid);
          if (muid != 0) {
            ValNodeAddInt (&pub, PUB_Muid, muid);
          }
        }
      }
      MonitorStrValue (mon, "Closing MedArch");
      Update ();
      MedArchFini ();
      MonitorFree (mon);
      ArrowCursor ();
      Update ();
    } else {
      ArrowCursor ();
      Message (MSG_ERROR, "Unable to connect to MedArch");
      MonitorFree (mon);
      Update ();
    }
  }
  return pub;
}


static Boolean HasMoreThanOneISOJTA (TitleMsgPtr titles)
{
  TitleMsgPtr tmp;
  ValNodePtr  ttl;
  Int4 num_found = 0;

  for (tmp = titles; tmp != NULL; tmp = tmp->next) {
    for (ttl = tmp->title; ttl != NULL; ttl = ttl->next) {
      if (ttl->choice == Cit_title_iso_jta) {
        if (num_found > 0) {
          return TRUE;
        }
        num_found++;
        break;
      }
    }
  }
  return FALSE;
}


static CharPtr GetExtendedDescription (TitleMsgPtr tmp) 
{
  CharPtr iso_jta = NULL, name = NULL, issn = NULL, extended = NULL;
  ValNodePtr ttl;

  if (tmp == NULL) {
    return NULL;
  }
  
  for (ttl = tmp->title; ttl != NULL; ttl = ttl->next) {
    if (ttl->choice == Cit_title_iso_jta) {
      iso_jta = ttl->data.ptrvalue;
    } else if (ttl->choice == Cit_title_name) {
      name = ttl->data.ptrvalue;
    } else if (ttl->choice == Cit_title_issn) {
      issn = ttl->data.ptrvalue;
    }
  }
  if (iso_jta == NULL) {
    return NULL;
  }
  if (name == NULL && issn == NULL) {
    extended = StringSave (iso_jta);
  } else if (name == NULL) {
    extended = (CharPtr) MemNew (sizeof (Char) * (StringLen (iso_jta) + StringLen (issn) + 5));
    sprintf (extended, "%s||(%s)", iso_jta, issn);
  } else if (issn == NULL) {
    extended = (CharPtr) MemNew (sizeof (Char) * (StringLen (iso_jta) + StringLen (name) + 5));
    sprintf (extended, "%s||(%s)", iso_jta, name);
  } else {
    extended = (CharPtr) MemNew (sizeof (Char) * (StringLen (iso_jta) + StringLen (name) + StringLen (issn) + 6));
    sprintf (extended, "%s||(%s:%s)", iso_jta, name, issn);
  }
  return extended;
}

static Boolean LookupJournalFuncNew (CharPtr title, size_t maxsize, Int1Ptr jtaType, ValNodePtr PNTR all_titlesP)

{
  ValNodePtr       first = NULL;
  ValNodePtr       last = NULL;
  MlaBackPtr       mbp;
  MlaRequestPtr    mrp;
  CharPtr          str, end;
  TitleMsgListPtr  tlp;
  TitleMsgPtr      tmp;
  ValNodePtr       ttl;
  ValNodePtr       vnp;

  if (all_titlesP != NULL) {
    *all_titlesP = NULL;
  }
  if (jtaType != NULL) {
    *jtaType = 0;
  }
  if (StringHasNoText (title)) return FALSE;

#ifdef OS_UNIX
  if (! log_mla_set) {
    str = (CharPtr) getenv ("LOG_MLA_ASN");
    if (StringDoesHaveText (str)) {
      if (StringICmp (str, "TRUE") == 0) {
        log_mla_asn = TRUE;
      }
    }
    log_mla_set = TRUE;
  }
#endif

  WatchCursor ();
  Update ();

  mrp = Mla2CreateJournalTitleRequest (title);
  if (mrp != NULL) {
  if (log_mla_asn) {
    LaunchAsnTextViewer ((Pointer) mrp, (AsnWriteFunc) MlaRequestAsnWrite, "lookup journal request");
  }
    mbp = Mla2SynchronousQuery (mrp);
    mrp = MlaRequestFree (mrp);
    if (mbp != NULL) {
      if (log_mla_asn) {
        LaunchAsnTextViewer ((Pointer) mbp, (AsnWriteFunc) MlaBackAsnWrite, "lookup journal result");
      }
      tlp = Mla2ExtractJournalTitleReply (mbp);
      if (tlp != NULL) {
        if (HasMoreThanOneISOJTA (tlp->titles)) {
          for (tmp = tlp->titles; tmp != NULL; tmp = tmp->next) {
            str = GetExtendedDescription (tmp);
            if (str != NULL) {
              ValNodeAddPointer (&first, 0, str);
            }
          }
        } else {
          for (tmp = tlp->titles; tmp != NULL; tmp = tmp->next) {
            for (ttl = tmp->title; ttl != NULL; ttl = ttl->next) {
              if (ttl->choice != Cit_title_iso_jta) continue;
              str = (CharPtr) ttl->data.ptrvalue;
              if (StringHasNoText (str)) continue;
              vnp = ValNodeCopyStr (&last, ttl->choice, str);
              if (first == NULL) {
                first = vnp;
              }
              last = vnp;
            }
          }
        }
        if (first != NULL) {
          str = (CharPtr) first->data.ptrvalue;
          StringNCpy_0 (title, str, maxsize);
          end = StringSearch (title, "||");
          if (end != NULL) {
            *end = 0;
          }
          if (jtaType != NULL) {
            *jtaType = Cit_title_iso_jta;
          }
        }
        if (all_titlesP == NULL) {
          first = ValNodeFreeData (first);
        } else {
          *all_titlesP = first;
        }
        tlp = TitleMsgListFree (tlp);
      }
      mbp = MlaBackFree (mbp);
    }
  }

  ArrowCursor ();
  Update ();

  return TRUE;
}

static Boolean LookupJournalFunc (CharPtr title, size_t maxsize, Int1Ptr jtaType, ValNodePtr PNTR all_titlesP)

{
  ValNodePtr  first = NULL;
  Int4        i;
  ValNodePtr  last = NULL;
  MonitorPtr  mon;
  Int4        num;
  Char        str [256];
  CharPtr     titles [20];
  Int1        types [20];
  ValNodePtr  vnp;

  WatchCursor ();
  mon = MonitorStrNewEx ("Lookup Journal", 40, FALSE);
  MonitorStrValue (mon, "Connecting to MedArch");
  Update ();
  if (MedArchInit ()) {
    StringNCpy_0 (str, title, sizeof (str));
    if (str [0] != '\0') {
      MonitorStrValue (mon, "Performing Lookup");
      Update ();
      num = MedArchGetTitles (titles, types, str, (Int1) Cit_title_jta,
                              Cit_title_iso_jta, 20);
      if (num > 0 && all_titlesP != NULL) {
        for (i = 0; i < num; i++) {
          vnp = ValNodeCopyStr (&last, types [i], titles [i]);
          if (first == NULL) {
            first = vnp;
          }
          last = vnp;
        }
        *all_titlesP = first;
      }
      if (num > 0 && types [0] == Cit_title_iso_jta) {
        StringNCpy_0 (title, titles [0], maxsize);
        if (jtaType != NULL) {
          *jtaType = types [0];
        }
        for (i = 0; i < num; i++) {
          MemFree (titles [i]);
        }
      } else {
        /*
        Message (MSG_OK, "Unable to match journal");
        */
      }
    }
    MonitorStrValue (mon, "Closing MedArch");
    Update ();
    MedArchFini ();
    MonitorFree (mon);
    ArrowCursor ();
    Update ();
    return TRUE;
  } else {
    ArrowCursor ();
    Message (MSG_ERROR, "Unable to connect to MedArch");
    MonitorFree (mon);
    Update ();
  }
  return FALSE;
}
/*#endif*/

/*#ifdef USE_TAXON*/
static Boolean LookupTaxonomyFunc (Uint2 entityID)

{
  SeqEntryPtr  sep;

  if (entityID < 1) return FALSE;
  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep == NULL) return FALSE;
  if (! leaveAsOldAsn) {
    MySeqEntryToAsn3 (sep, TRUE, FALSE, TRUE);
  }
  return TRUE;
}
/*#endif*/

extern void QuitProc (void);
extern void QuitProc (void)

{
  Boolean        dirty;
  Uint4          j;
  Uint4          num;
  ObjMgrPtr      omp;
  ObjMgrDataPtr  PNTR omdpp;
  ObjMgrDataPtr  tmp;
#ifdef USE_SMARTNET
  SMUserDataPtr  sm_usr_data = NULL;
  OMUserDataPtr  omudp;
#endif

#ifdef WIN_MAC
  if (initialFormsActive) {
    if (Message (MSG_YN,
        "This will end your session without saving any data.\nAre you sure you want to exit?") == ANS_YES) {
      QuitProgram ();
    }
    return;
  }
#endif

  if (smartnetMode) {
#ifdef USE_SMARTNET
      
      omp = ObjMgrGet ();
      num = omp->HighestEntityID;
      for(j = 1; j <= omp->HighestEntityID; j++) {
          if((omudp = ObjMgrGetUserData(j, 0,0, SMART_KEY)) != NULL) {
            if((sm_usr_data = 
                (SMUserDataPtr) omudp->userdata.ptrvalue) != NULL &&
               sm_usr_data->fd != 0) {
              sm_usr_data->header->status = SMStatClosed;
              SMSendMsgToClient(sm_usr_data);
            }              
          }
      }
#endif
      QuitProgram ();
      return;
  }

  dirty = FALSE;
  omp = ObjMgrGet ();
  num = omp->currobj;
  for (j = 0, omdpp = omp->datalist; j < num && omdpp != NULL; j++, omdpp++) {
    tmp = *omdpp;
    if (tmp->parentptr == NULL) {
      if (tmp->dirty) {
        dirty = TRUE;
      }
    }
  }
  if (dirty) {
    if (Message (MSG_YN,
        "Some data have not been saved.\nAre you sure you want to exit?") == ANS_NO) {
      return;
    }
  } else if (subtoolMode || smartnetMode) {
    if (Message (MSG_YN,
        "Are you sure you want to exit?") == ANS_NO) {
      return;
    }
  }
#ifndef WIN_MAC
  if (subtoolMode || smartnetMode || backupMode) {
    subtoolRecordDirty = FALSE;
    FileRemove (SEQUIN_EDIT_TEMP_FILE);
    /* FileRemove (SEQUIN_EDIT_PREV_FILE); */
    FileRemove (SEQUIN_EDIT_BACK_FILE);
    FileRemove (SEQUIN_EDIT_ARCH_FILE);
    FileRename (SEQUIN_EDIT_PREV_FILE, SEQUIN_EDIT_ARCH_FILE);
    if (backupMode) {
      FileRemove (SEQUIN_EDIT_ARCH_FILE);
    }
  }
#endif
  QuitProgram ();
}

static void CloseProc (BaseFormPtr bfp)

{
  BioseqViewFormPtr  bvfp;
  Uint2              entityID;
  Boolean            freeEditors = FALSE;
  Int4               num_dirty;
  Uint4              j;
  Uint4              num;
  Boolean            numview;
  ObjMgrPtr          omp;
  ObjMgrDataPtr      PNTR omdpp;
  OMUserDataPtr      omudp;
  ObjMgrDataPtr      tmp;
  Uint2              userkey;
  ValNodePtr         vnp;
#ifdef USE_SMARTNET
  SMUserDataPtr     sm_usr_data = NULL;
#endif

#ifdef WIN_MAC
  if (initialFormsActive) return;
#endif
  if (bfp != NULL) {
    omp = ObjMgrGet ();
    num = omp->currobj;
    for (j = 0, omdpp = omp->datalist; j < num && omdpp != NULL; j++, omdpp++) {
      tmp = *omdpp;
      if (tmp->parentptr == NULL && tmp->EntityID == bfp->input_entityID) {
        num_dirty = 0;
        numview = 0;
        for (omudp = tmp->userdata; omudp != NULL; omudp = omudp->next) {
          if (omudp->proctype == OMPROC_VIEW) {
            numview++;
            if (tmp->dirty) {
              num_dirty++;
            }
          }
        }
        if (num_dirty == 1 && numview < 2) {
          if (Message (MSG_OKC,
              "Closing the window will lose unsaved data.") != ANS_OK) {
            return;
          }
        }
        if (numview < 2) {
          freeEditors = TRUE;
        }
      }
    }
    if (freeEditors) {
      bvfp = (BioseqViewFormPtr) bfp;
      if (bvfp->input_entityID > 0 && bvfp->userkey > 0) {
        /* need to unlock far components before freeing data */
        if (bvfp->bvd.bsplist != NULL) {
          bvfp->bvd.bsplist = UnlockFarComponents (bvfp->bvd.bsplist);
        }
        userkey = bvfp->userkey;
        bvfp->userkey = 0;
        /* ObjMgrFreeUserData (bvfp->input_entityID, bvfp->procid, bvfp->proctype, userkey); */
        /* this may trigger another remove, hence bvfp->userkey first set to 0 */
        for (vnp = bvfp->bvd.entityList; vnp != NULL; vnp = vnp->next) {
          if (bvfp->input_entityID != (Uint2) vnp->data.intvalue) {
            ObjMgrFreeUserData ((Uint2) vnp->data.intvalue, bvfp->procid, bvfp->proctype, userkey);
          }
        }
        bvfp->userkey = userkey;
      }
      if (bvfp->bvd.entityList != NULL) {
        bvfp->bvd.entityList = ValNodeFree (bvfp->bvd.entityList);
      }
    }
    numview = 0;
    for (j = 0, omdpp = omp->datalist; j < num && omdpp != NULL; j++, omdpp++) {
      tmp = *omdpp;

      if (tmp->parentptr == NULL) {
          for (omudp = tmp->userdata; omudp != NULL; omudp = omudp->next) {

              if (omudp->proctype == OMPROC_VIEW) {
                  numview++;
              }
          }
      }
    }

    entityID = bfp->input_entityID;
    if(!smartnetMode) {
        /* RemoveSeqEntryViewer (bfp->form); */ /* can go back to Remove */
        Remove (bfp->form);
        if (freeEditors) {
          ObjMgrSendMsg (OM_MSG_DEL, entityID, 0, 0);
        }
    }

    if (numview <= 1) {
        if (subtoolMode || stdinMode || binseqentryMode) {
            Message (MSG_OK, "No more viewers, quitting program.");
            QuitProgram ();
            return;
        } else if (smartnetMode) {
#ifdef USE_SMARTNET
          if((omudp = 
              ObjMgrGetUserData(bfp->input_entityID, 0, 0, SMART_KEY)) != NULL) {
              if((sm_usr_data = 
                  (SMUserDataPtr) omudp->userdata.ptrvalue) != NULL && 
                 sm_usr_data->fd != (int) NULL) {
                  sm_usr_data->header->status = SMStatClosed;
                  SMSendMsgToClient(sm_usr_data); 
              }

              RemoveSeqEntryViewer (bfp->form);
              ObjMgrFreeUserData(entityID, 0, 0, SMART_KEY); 
              ObjMgrFreeUserData (entityID, 0, 0, DUMB_KEY);
              if (freeEditors) {
                ObjMgrSendMsg (OM_MSG_DEL, entityID, 0, 0);
                if (DeleteRemainingViews(entityID)) {
                  VSeqMgrShow();
                }
              }
              
          }
          return;
#endif
        }

        WatchCursor ();
        Update ();
        if (! workbenchMode) {
            if (termListForm != NULL || docSumForm != NULL) {
            } else {
                Show (startupForm);
                Select (startupForm);
                SendHelpScrollMessage (helpForm, "Initial Entry", NULL);
            }
            ArrowCursor ();
            Update ();
        }
    } else { /* numview > 1 */
        if(smartnetMode) {
            entityID = bfp->input_entityID;
            RemoveSeqEntryViewer (bfp->form);
#ifdef USE_SMARTNET 
            /* sssddd   ObjMgrFreeUserData(entityID, 0, 0, SMART_KEY); */
#endif 
            if (freeEditors) {
              ObjMgrSendMsg (OM_MSG_DEL, entityID, 0, 0);
            }
        }

    }
  }
}

static void EditSubmitBlock (BaseFormPtr bfp)

{
  Int2           handled;
  ObjMgrDataPtr  omdp;

  if (bfp != NULL && bfp->input_entityID != 0) {
    omdp = ObjMgrGetData (bfp->input_entityID);
    if (omdp != NULL) {
      if (omdp->datatype == OBJ_SEQSUB) {
        WatchCursor ();
        handled = GatherProcLaunch (OMPROC_EDIT, FALSE, bfp->input_entityID, 1,
                                    OBJ_SUBMIT_BLOCK, 0, 0, OBJ_SUBMIT_BLOCK, 0);
        ArrowCursor ();
        if (handled != OM_MSG_RET_DONE || handled == OM_MSG_RET_NOPROC) {
          Message (MSG_FATAL, "Unable to launch editor.");
        }
      } else {
        Message (MSG_OK, "Record has no submit block.  This may be added when record is first read.");
      }
    }
  }
}

/*
static void DisplayFontChangeProc (IteM i)

{
  BaseFormPtr  bfp;
  FonT         fnt;
  FontSpec     fs;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp != NULL) {
    MemSet ((Pointer) (&fs), 0, sizeof (FontSpec));
    fnt = programFont;
    if (seqviewprocs.displayFont != NULL) {
      fnt = seqviewprocs.displayFont;
    }
    GetFontSpec (fnt, &fs);
    if (ChooseFont (&fs, CFF_READ_FSP | CFF_MONOSPACE, NULL)) {
      seqviewprocs.displayFont = CreateFont (&fs);
      SendMessageToForm (bfp->form, VIB_MSG_REDRAW);
    }
  }
}
*/

static void DuplicateViewProc (IteM i)

{
  BaseFormPtr  bfp;
  Int2         handled;
  Uint4        itemID;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  if (bfp->input_itemtype == OBJ_BIOSEQ) {
    WatchCursor ();
    itemID = bfp->input_itemID;
    if (itemID == 0) {
      itemID = 1;
    }
    seqviewprocs.forceSeparateViewer = TRUE;
    SeqEntrySetScope (NULL);
    handled = GatherProcLaunch (OMPROC_VIEW, FALSE, bfp->input_entityID, itemID,
                                OBJ_BIOSEQ, 0, 0, OBJ_BIOSEQ, 0);
    ArrowCursor ();
    if (handled != OM_MSG_RET_DONE || handled == OM_MSG_RET_NOPROC) {
      Message (MSG_ERROR, "Unable to launch additional viewer.");
    }
  } else if (bfp->input_itemtype == OBJ_MEDLINE_ENTRY) {
    WatchCursor ();
    itemID = bfp->input_itemID;
    if (itemID == 0) {
      itemID = 1;
    }
    handled = GatherProcLaunch (OMPROC_VIEW, FALSE, bfp->input_entityID, itemID,
                                OBJ_MEDLINE_ENTRY, 0, 0, OBJ_MEDLINE_ENTRY, 0);
    ArrowCursor ();
    if (handled != OM_MSG_RET_DONE || handled == OM_MSG_RET_NOPROC) {
      Message (MSG_ERROR, "Unable to launch additional viewer.");
    }
  }
}


static void RestoreSeqEntryProc (IteM i)

{
  BaseFormPtr  bfp;
  Char         path [PATH_MAX];
  Uint2        new_entityID;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp != NULL && bfp->input_itemtype == OBJ_BIOSEQ) {
    if (GetInputFileName (path, sizeof (path), "", "TEXT")) {
      new_entityID = RestoreEntityIDFromFile (path, bfp->input_entityID);
      if (new_entityID != 0) {
        bfp->input_entityID = new_entityID;
        ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
        ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
      }
    }
  }
}

void AddAboutAndHelpMenuItems (MenU m)

{
  CommandItem (m, "About Sequin...", AboutProc);
  CommandItem (m, "Help...", DisplayHelpFormProc);
  SeparatorItem (m);
}

static void FindOrf (IteM i)

{
  BaseFormPtr  bfp;
  BioseqPtr    bsp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  bsp =  GetBioseqGivenIDs (bfp->input_entityID, bfp->input_itemID, bfp->input_itemtype);
  if (bsp == NULL) return;
  LaunchOrfViewer (bsp, bfp->input_entityID, bfp->input_itemID, FALSE);
}

static void FindAlu (IteM i)

{
  MsgAnswer    ans;
  BaseFormPtr  bfp;
  BioseqPtr    bsp;
  Boolean      got_some;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  bsp =  GetBioseqGivenIDs (bfp->input_entityID, bfp->input_itemID, bfp->input_itemtype);
  if (bsp == NULL) return;
  if (bsp->length >= 100) {
    ans = Message (MSG_OKC, "This calculation may take a while");
    if (ans == ANS_CANCEL) return;
  }
  WatchCursor ();
  got_some = FindHumanRepeats (bsp, TRUE);
  ArrowCursor ();
  if (got_some) {
    Message (MSG_OK, "Repeat regions were found.");
    ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
    ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
  } else {
    Message (MSG_OK, "No Repeat regions found.");
  }
}

typedef struct changetargetform {
  FORM_MESSAGE_BLOCK
  TexT               seqid;
  BaseFormPtr        base;
} ChangeTargetForm, PNTR ChangeTargetFormPtr;

static void AcceptChangeTargetProc (ButtoN b)

{
  ChangeTargetFormPtr  cfp;
  Char                 str [128];

  cfp = (ChangeTargetFormPtr) GetObjectExtra (b);
  if (cfp == NULL) return;
  GetTitle (cfp->seqid, str, sizeof (str));
  SetBioseqViewTarget (cfp->base, str);
}

static void ClearTargetBtn (ButtoN b)

{
  ChangeTargetFormPtr  cfp;

  cfp = (ChangeTargetFormPtr) GetObjectExtra (b);
  if (cfp == NULL) return;
  SetTitle (cfp->seqid, "");
}

static void ChangeTargetMessageProc (ForM f, Int2 mssg)

{
  StdEditorProcsPtr  sepp;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    if (sepp->handleMessages != NULL) {
      sepp->handleMessages (f, mssg);
    }
  }
}

extern void ChangeTargetBaseForm (BaseFormPtr bfp)
{
  ButtoN               b;
  GrouP                c;
  ChangeTargetFormPtr  cfp;
  GrouP                g;
  StdEditorProcsPtr    sepp;
  WindoW               w;

  if (bfp == NULL) return;

  cfp = (ChangeTargetFormPtr) MemNew (sizeof (ChangeTargetForm));
  if (cfp == NULL) return;
  w = FixedWindow (-50, -33, -10, -10, "Select Target by SeqID", StdCloseWindowProc);
  SetObjectExtra (w, cfp, StdCleanupFormProc);
  cfp->form = (ForM) w;
  cfp->formmessage = ChangeTargetMessageProc;

  g = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (g, 10, 10);

  cfp->seqid = DialogText (g, "", 12, NULL);
  cfp->base = bfp;
  c = HiddenGroup (g, 3, 0, NULL);
  SetGroupSpacing (c, 10, 10);
  b = DefaultButton (c, "Accept", AcceptChangeTargetProc);
  SetObjectExtra (b, cfp, NULL);
  b = PushButton (c, "Clear", ClearTargetBtn);
  SetObjectExtra (b, cfp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);

  AlignObjects (ALIGN_CENTER, (HANDLE) cfp->seqid, (HANDLE) c, NULL);
  RealizeWindow (w);
  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    SetActivate (w, sepp->activateForm);
  }
  Show (w);
  Select (w);
}

static void DoChangeTarget (IteM i)

{
  BaseFormPtr          bfp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  ChangeTargetBaseForm (bfp);
}

static void ConfigFormMessage (ForM f, Int2 mssg)

{
  BaseFormPtr        bfp;
  StdEditorProcsPtr  sepp;

  bfp = (BaseFormPtr) GetObjectExtra (f);
  if (bfp != NULL) {
    switch (mssg) {
      case VIB_MSG_CUT :
        StdCutTextProc (NULL);
        break;
      case VIB_MSG_COPY :
        StdCopyTextProc (NULL);
        break;
      case VIB_MSG_PASTE :
        StdPasteTextProc (NULL);
        break;
      case VIB_MSG_DELETE :
        StdDeleteTextProc (NULL);
        break;
      default :
        sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
        if (sepp != NULL && sepp->handleMessages != NULL) {
          sepp->handleMessages (f, mssg);
        }
        break;
    }
  }
}

static void ConfigAccepted (void)

{
  if (WriteSequinAppParam ("SETTINGS", "PUBLICNETWORKSEQUIN", "TRUE")) {
    Message (MSG_OK, "Setting will take affect when you restart Sequin");
  }
}

static void ConfigCancelled (void)

{
  Message (MSG_OK, "No changes to the network configuration have been made");
}

static void ConfigTurnedOff (void)

{
  if (WriteSequinAppParam ("SETTINGS", "PUBLICNETWORKSEQUIN", "FALSE")) {
    Message (MSG_OK, "Setting will take affect when you restart Sequin");
  }
}

typedef struct prefsform {
  FORM_MESSAGE_BLOCK
  DialoG             prefs;
} PrefsForm, PNTR PrefsFormPtr;

/*
static void AcceptPrefsProc (ButtoN b)

{
  EntrezPrefsPtr  epp;
  PrefsFormPtr    pfp;

  pfp = (PrefsFormPtr) GetObjectExtra (b);
  if (pfp == NULL) return;
  Hide (pfp->form);
  epp = (EntrezPrefsPtr) DialogToPointer (pfp->prefs);
  if (epp != NULL) {
    entrezglobals.persistDefault = epp->persistDefault;
    entrezglobals.alignDefault = epp->alignDefault;
    entrezglobals.lookupDirect = epp->lookupDirect;
    entrezglobals.showAsn = epp->showAsn;
    ReplaceString (&entrezglobals.initDatabase, epp->initDatabase);
    ReplaceString (&entrezglobals.initField, epp->initField);
    ReplaceString (&entrezglobals.initMode, epp->initMode);
    ReplaceString (&medviewprocs.initMedLabel, epp->initMedLabel);
    ReplaceString (&seqviewprocs.initNucLabel, epp->initNucLabel);
    ReplaceString (&seqviewprocs.initProtLabel, epp->initProtLabel);
    ReplaceString (&seqviewprocs.initGenomeLabel, epp->initGenomeLabel);
    SetSequinAppParamTF ("PREFERENCES", "PARENTSPERSIST", entrezglobals.persistDefault);
    SetSequinAppParamTF ("PREFERENCES", "ALIGNCHECKED", entrezglobals.alignDefault);
    SetSequinAppParamTF ("PREFERENCES", "LOOKUPDIRECT", entrezglobals.lookupDirect);
    SetSequinAppParamTF ("PREFERENCES", "SHOWASNPAGE", entrezglobals.showAsn);
    SetSequinAppParam ("SETTINGS", "DATABASE", entrezglobals.initDatabase);
    SetSequinAppParam ("SETTINGS", "FIELD", entrezglobals.initField);
    SetSequinAppParam ("SETTINGS", "MODE", entrezglobals.initMode);
    SetSequinAppParam ("SETTINGS", "MEDPAGE", medviewprocs.initMedLabel);
    SetSequinAppParam ("SETTINGS", "NUCPAGE", seqviewprocs.initNucLabel);
    SetSequinAppParam ("SETTINGS", "PRTPAGE", seqviewprocs.initProtLabel);
    SetSequinAppParam ("SETTINGS", "GENMPAGE", seqviewprocs.initGenomeLabel);
  }
  EntrezPrefsFree (epp);
  Remove (pfp->form);
}
*/

static void DefaultMessageProc (ForM f, Int2 mssg)

{
  StdEditorProcsPtr  sepp;

  sepp = (StdEditorProcsPtr) GetAppProperty ("StdEditorForm");
  if (sepp != NULL) {
    if (sepp->handleMessages != NULL) {
      sepp->handleMessages (f, mssg);
    }
  }
}

static void PreferencesProc (IteM i)

{
  /*
  ButtoN          b;
  GrouP           c;
  EntrezPrefsPtr  epp;
  GrouP           g;
  PrefsFormPtr    pfp;
  WindoW          w;

  pfp = (PrefsFormPtr) MemNew (sizeof (PrefsForm));
  if (pfp == NULL) return;
  if (! EntrezIsInited ()) {
    EntrezBioseqFetchEnable ("Sequin", TRUE);
    SequinEntrezInit ("Sequin", FALSE, NULL);
  }
  w = FixedWindow (-50, -33, -10, -10, "Preferences", StdCloseWindowProc);
  SetObjectExtra (w, pfp, StdCleanupFormProc);
  pfp->form = (ForM) w;
  pfp->formmessage = DefaultMessageProc;
  g = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (g, 10, 10);
  pfp->prefs = CreateEntrezPrefsDialog (g, NULL);
  c = HiddenGroup (g, 2, 0, NULL);
  b = DefaultButton (c, "Accept", AcceptPrefsProc);
  SetObjectExtra (b, pfp, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);
  AlignObjects (ALIGN_CENTER, (HANDLE) pfp->prefs, (HANDLE) c, NULL);
  RealizeWindow (w);
  epp = EntrezPrefsNew ();
  if (epp != NULL) {
    epp->persistDefault = entrezglobals.persistDefault;
    epp->alignDefault = entrezglobals.alignDefault;
    epp->lookupDirect = entrezglobals.lookupDirect;
    epp->showAsn = entrezglobals.showAsn;
    epp->initDatabase = StringSaveNoNull (entrezglobals.initDatabase);
    epp->initField = StringSaveNoNull (entrezglobals.initField);
    epp->initMode = StringSaveNoNull (entrezglobals.initMode);
    epp->initMedLabel = StringSaveNoNull (medviewprocs.initMedLabel);
    epp->initNucLabel = StringSaveNoNull (seqviewprocs.initNucLabel);
    epp->initProtLabel = StringSaveNoNull (seqviewprocs.initProtLabel);
    epp->initGenomeLabel = StringSaveNoNull (seqviewprocs.initGenomeLabel);
  }
  PointerToDialog (pfp->prefs, (Pointer) epp);
  EntrezPrefsFree (epp);
  Show (w);
  Select (w);
  */
}

void NetConfigureProc (IteM i)

{
  Boolean  netCurrentlyOn = FALSE;
  Char     str [32];

  if (GetSequinAppParam ("SETTINGS", "PUBLICNETWORKSEQUIN", NULL, str, sizeof (str))) {
    if (StringICmp (str, "TRUE") == 0) {
      netCurrentlyOn = TRUE;
    }
  }
  if (useEntrez) {
    netCurrentlyOn = TRUE;
  }
  SendHelpScrollMessage (helpForm, "Misc Menu", "Net Configure");
  ShowNetConfigForm (ConfigFormActivated, ConfigFormMessage,
                     ConfigAccepted, ConfigCancelled,
                     ConfigTurnedOff, netCurrentlyOn);
}


extern void SqnReadAlignView (BaseFormPtr bfp, BioseqPtr target_bsp, SeqEntryPtr source_sep, Boolean do_update)

{
  BioseqPtr  nbsp;

  if (target_bsp == NULL || source_sep == NULL) return;
  nbsp = FindNucBioseq (source_sep);
  if (nbsp == NULL) return;

  if (do_update) {
    UpdateSeqAfterDownload (bfp, target_bsp, nbsp);
  } else {
    ExtendSeqAfterDownload (bfp, target_bsp, nbsp);
  }
}

extern void NewFeaturePropagate (
  IteM i
);

static void ExtendSeqWithFASTA (IteM i)

{
  NewExtendSequence (i);
}

static void ExtendSeqWithRec (IteM i)

{
  NewExtendSequence (i);
}

NLM_EXTERN Int4 GenomeFromLocName (CharPtr loc_name);

/*******COLOMBE ***********/


static void CommonAddSeq (IteM i, Int2 type)

{
  Nlm_QualNameAssocPtr qp;
  BaseFormPtr        bfp;
  Uint1              biomol;
  BioSourcePtr       biop = NULL;
  BioseqPtr          bsp;
  BioseqSetPtr       bssp;
  Pointer            dataptr;
  Uint2              datatype;
  FILE               *fp;
  MolInfoPtr         mip;
  MolInfoPtr         molinf = NULL;
  OrgRefPtr          orp;
  SeqEntryPtr        nuc;
  Char               path [PATH_MAX];
  CharPtr            ptr;
  SeqEntryPtr        sep;
  BioSourcePtr       src = NULL;
  Char               str [128];
  CharPtr            tax = NULL;
  CharPtr            title = NULL;
  SeqEntryPtr        top;
  ValNodePtr         vnp;
  Int4               tmp;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  top = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (top == NULL) return;
  if (type == 1 || type == 2) {
    if (! GetInputFileName (path, sizeof (path),"","TEXT")) return;
    fp = FileOpen (path, "r");
    if (fp == NULL) return;
    if (IS_Bioseq_set (top)) {
      bssp = (BioseqSetPtr) top->data.ptrvalue;
      if (bssp != NULL && bssp->_class != BioseqseqSet_class_phy_set) {
        sep = bssp->seq_set;
        vnp = SeqEntryGetSeqDescr (top, Seq_descr_source, NULL);
        if (vnp != NULL) {
          src = (BioSourcePtr) vnp->data.ptrvalue;
          if (src != NULL) {
            orp = biop->org;
            if (orp != NULL) {
              tax = orp->taxname;
            }
          }
        }
      }
    }
    nuc = FindNucSeqEntry (top);
    vnp = SeqEntryGetSeqDescr (nuc, Seq_descr_molinfo, NULL);
    if (vnp != NULL) {
      molinf = (MolInfoPtr) vnp->data.ptrvalue;
    }
    while ((dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, NULL, FALSE, FALSE, TRUE, FALSE)) != NULL) {
      if (datatype == OBJ_BIOSEQ || datatype == OBJ_BIOSEQSET) {
        bsp = NULL;
        bssp = NULL;
        sep = SeqMgrGetSeqEntryForData (dataptr);
        if (sep == NULL) {
          sep = SeqEntryNew ();
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
        if (sep != NULL) {
          AddSeqEntryToSeqEntry (top, sep, TRUE);
          title = SeqEntryGetTitle (sep);
          vnp = SeqEntryGetSeqDescr (sep, Seq_descr_source, NULL);
          if (vnp == NULL || title != NULL) {
            ptr = StringISearch (title, "[org=");
            if (ptr != NULL) {
              StringNCpy_0 (str, ptr + 5, sizeof (str));
              ptr = StringChr (str, ']');
            } else {
              ptr = StringISearch (title, "[organism=");
              if (ptr != NULL) {
                StringNCpy_0 (str, ptr + 10, sizeof (str));
                ptr = StringChr (str, ']');
              }
            }
            if (ptr != NULL) {
              *ptr = '\0';
              biop = BioSourceNew ();
              if (biop != NULL) {
                orp = OrgRefNew ();
                biop->org = orp;
                if (orp != NULL) {
                  SetTaxNameAndRemoveTaxRef (orp, StringSave (str));
                }
                vnp = CreateNewDescriptor (sep, Seq_descr_source);
                if (vnp != NULL) {
                  vnp->data.ptrvalue = (Pointer) biop;
                }
              }
            }
          }
          if (vnp == NULL && tax != NULL) {
            biop = BioSourceNew ();
            if (biop != NULL) {
              orp = OrgRefNew ();
              biop->org = orp;
              if (orp != NULL) {
                SetTaxNameAndRemoveTaxRef (orp, tax);
              }
              vnp = CreateNewDescriptor (sep, Seq_descr_source);
              if (vnp != NULL) {
                vnp->data.ptrvalue = (Pointer) biop;
              }
            }
          }
          vnp = SeqEntryGetSeqDescr (sep, Seq_descr_source, NULL);
          if (vnp != NULL) {
            biop = (BioSourcePtr) vnp->data.ptrvalue;
          }
          if (biop != NULL && title != NULL) {
            for (qp = current_orgmod_subtype_alist; qp->name != NULL; qp++) {
              MakeSearchStringFromAlist (str, qp->name);
              AddToOrgMod (biop, title, str, qp->value);
              ExciseString (title, str, "]");
            }
            for (qp = current_subsource_subtype_alist; qp->name != NULL; qp++) {
              MakeSearchStringFromAlist (str, qp->name);
              AddToSubSource (biop, title, str, qp->value);
              ExciseString (title, str, "]");
            }
            AddToOrgMod (biop, title, "[note=", 255);
            ExciseString (title, "[note=", "]");
            AddToSubSource (biop, title, "[subsource=", 255);
            ExciseString (title, "[subsource=", "]");
            ExciseString (title, "[org=", "]");
            ExciseString (title, "[organism=", "]");
            if (bsp != NULL) {
              ptr = StringISearch (title, "[molecule=");
              if (ptr != NULL) {
                StringNCpy_0 (str, ptr + 10, sizeof (str));
                ptr = StringChr (str, ']');
                if (ptr != NULL) {
                  *ptr = '\0';
                  if (StringCmp (str, "dna") == 0) {
                    bsp->mol = Seq_mol_dna;
                  } else if (StringCmp (str, "rna") == 0) {
                    bsp->mol = Seq_mol_rna;
                  }
                }
              }
            }
            ptr = StringISearch (title, "[location=");
            if (ptr != NULL) {
              StringNCpy_0 (str, ptr + 10, sizeof (str));
              ptr = StringChr (str, ']');
              if (ptr != NULL) {
                *ptr = '\0';
                if (StringICmp (str, "Mitochondrial") == 0) { /* alternative spelling */
                  biop->genome = 5;
                }
                tmp = GenomeFromLocName (str);
                if (tmp > -1) {
                  biop->genome = tmp;
                }
              }
            }
            ptr = StringISearch (title, "[moltype=");
            if (ptr != NULL) {
              biomol = 0;
              StringNCpy_0 (str, ptr + 8, sizeof (str));
              ptr = StringChr (str, ']');
              if (ptr != NULL) {
                *ptr = '\0';
                if (StringICmp (str, "genomic") == 0) {
                  biomol = MOLECULE_TYPE_GENOMIC;
                } else if (StringICmp (str, "mRNA") == 0) {
                  biomol = MOLECULE_TYPE_MRNA;
                }
                if (biomol != 0) {
                  mip = MolInfoNew ();
                  if (mip != NULL) {
                    mip->biomol = biomol;
                    vnp = CreateNewDescriptor (sep, Seq_descr_molinfo);
                    if (vnp != NULL) {
                      vnp->data.ptrvalue = (Pointer) mip;
                    }
                  }
                }
              }
            }
          }
          ExciseString (title, "[molecule=", "]");
          ExciseString (title, "[moltype=", "]");
          ExciseString (title, "[location=", "]");
          TrimSpacesAroundString (title);
          if (title != NULL && StringHasNoText (title)) {
            vnp = NULL;
            if (IS_Bioseq (sep)) {
              bsp = (BioseqPtr) sep->data.ptrvalue;
              vnp = ValNodeExtract (&(bsp->descr), Seq_descr_title);
            } else if (IS_Bioseq_set (sep)) {
              bssp = (BioseqSetPtr) sep->data.ptrvalue;
              vnp = ValNodeExtract (&(bssp->descr), Seq_descr_title);
            }
            if (vnp != NULL && StringHasNoText ((CharPtr) vnp->data.ptrvalue)) {
              vnp = ValNodeFreeData (vnp);
            }
          }
          sep = FindNucSeqEntry (sep);
          vnp = SeqEntryGetSeqDescr (sep, Seq_descr_molinfo, NULL);
          if (vnp == NULL && molinf != NULL) {
            mip = MolInfoNew ();
            if (mip != NULL) {
              vnp = CreateNewDescriptor (sep, Seq_descr_molinfo);
              if (vnp != NULL) {
                vnp->data.ptrvalue = (Pointer) mip;
                mip->biomol = molinf->biomol;
                mip->tech = molinf->tech;
                mip->completeness = molinf->completeness;
              }
            }
          }
        }
      }
    }
    FileClose (fp);
    ObjMgrSetDirtyFlag (bfp->input_entityID, TRUE);
    ObjMgrSendMsg (OM_MSG_UPDATE, bfp->input_entityID, 0, 0);
    Update ();
  } else {
    Message (MSG_OK, "Not yet implemented");
  }
}

static void AddSeqWithFASTA (IteM i)

{
  CommonAddSeq (i, 1);
}

static void AddSeqWithRec (IteM i)

{
  CommonAddSeq (i, 2);
}

/*#ifdef ALLOW_DOWNLOAD*/
BioseqPtr  updateTargetBspKludge = NULL;

static void ExtendSeqWithAcc (IteM i)

{
  BaseFormPtr  bfp;
  BioseqPtr    bsp;
  SeqEntryPtr  sep;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;
  bsp = GetBioseqGivenIDs (bfp->input_entityID, bfp->input_itemID, bfp->input_itemtype);
  if (bsp == NULL) return;
  updateTargetBspKludge = bsp;
  CommonFetchFromNet (DownloadAndExtendProc, StdCancelButtonProc);
}
/*#endif*/


#ifndef WIN_MAC
extern void CreateNewLayoutMenu (MenU m, BaseFormPtr bp);

static void MedlineViewFormMenus (WindoW w)

{
  BaseFormPtr  bfp;
  IteM         i;
  MenU         m;

  bfp = (BaseFormPtr) GetObjectExtra (w);
  if (bfp != NULL) {
    m = PulldownMenu (w, "File");
    FormCommandItem (m, "Close", bfp, VIB_MSG_CLOSE);
    SeparatorItem (m);
    i = CommandItem (m, "Duplicate", DuplicateViewProc);
    SetObjectExtra (i, bfp, NULL);
    SeparatorItem (m);
    FormCommandItem (m, "Export...", bfp, VIB_MSG_EXPORT);
    SeparatorItem (m);
    /*
    FormCommandItem (m, "Save", bfp, VIB_MSG_SAVE);
    FormCommandItem (m, "Save As...", bfp, VIB_MSG_SAVE_AS);
    SeparatorItem (m);
    */
    FormCommandItem (m, "Print...", bfp, VIB_MSG_PRINT);

    m = PulldownMenu (w, "Edit");
    FormCommandItem (m, COPY_MENU_ITEM, bfp, VIB_MSG_COPY);
  }
}


static void BioseqViewFormMenus (WindoW w)

{
  BaseFormPtr    bfp;
  IteM           i;
  MenU           m;
  Int2           mssgadd;
  Int2           mssgalign;
  Int2           mssgdup;
  Int2           mssgfeatprop;
  Int2           mssgseq;
  Int2           mssgsub;
  Int2           mssgupd, mssgupd_idx;
  Int2           mssgext;
  ObjMgrDataPtr  omdp;
  MenU           sub;

  bfp = (BaseFormPtr) GetObjectExtra (w);
  if (bfp != NULL) {
    m = PulldownMenu (w, "File/ F");
/*#ifdef INTERNAL_NCBI_SEQUIN*/
    if (indexerVersion) {
      if (subtoolMode || stdinMode || binseqentryMode) {
        FormCommandItem (m, "Abort Session", bfp, VIB_MSG_QUIT);
        SeparatorItem (m);
        CommandItem (m, "Accept Changes", SubtoolDoneProc);
        SeparatorItem (m);
      } else if (smartnetMode) {
#ifdef USE_SMARTNET
        FormCommandItem (m, "Abort Session", bfp, VIB_MSG_CLOSE);
        SeparatorItem (m);
        i = CommandItem (m, "Accept Changes", SmartnetDoneProc);
        SetObjectExtra (i, bfp, NULL);
        SeparatorItem (m);
#endif
      }
    }
/*#endif*/
    AddAboutAndHelpMenuItems (m);
    i = CommandItem (m, "Open...", ReadNewAsnProc);
    SetObjectExtra (i, bfp, NULL);
    if (indexerVersion) {
      i = CommandItem (m, "FASTA Nucleotide Direct to Sequence Editor", FastaNucDirectToSeqEdProc);
      SetObjectExtra (i, bfp, NULL);
    }
    SeparatorItem (m);
    FormCommandItem (m, "Close", bfp, VIB_MSG_CLOSE);
    SeparatorItem (m);
    /* FormCommandItem (m, "Import Nucleotide FASTA...", bfp, VIB_MSG_IMPORT); */
    FormCommandItem (m, "Export...", bfp, VIB_MSG_EXPORT);
    SeparatorItem (m);
    i = CommandItem (m, "Duplicate View", DuplicateViewProc);
    SetObjectExtra (i, bfp, NULL);
    SeparatorItem (m);
    FormCommandItem (m, "Save", bfp, VIB_MSG_SAVE);
    FormCommandItem (m, "Save As...", bfp, VIB_MSG_SAVE_AS);
    i = CommandItem (m, "Save As Binary Seq-entry...", SaveBinSeqEntry);
    SetObjectExtra (i, bfp, NULL);
    SeparatorItem (m);
    i = CommandItem (m, "Restore...", RestoreSeqEntryProc);
    SetObjectExtra (i, bfp, NULL);
    if ((! subtoolMode) && (! stdinMode) &&
        (! binseqentryMode) && (! smartnetMode)) {
      SeparatorItem (m);
      i = CommandItem (m, "Prepare Submission...", PrepareSeqSubmitProc);
      SetObjectExtra (i, bfp, NULL);
      omdp = ObjMgrGetData (bfp->input_entityID);
      if (omdp != NULL && omdp->datatype != OBJ_SEQSUB) {
        Disable (i);
      }
      /*
      i = CommandItem (m, "Submit to NCBI", SubmitToNCBI);
      SetObjectExtra (i, bfp, NULL);
      if (omdp != NULL && omdp->datatype != OBJ_SEQSUB) {
        Disable (i);
      }
      */
    }
/*#ifdef INTERNAL_NCBI_SEQUIN*/
    if (indexerVersion) {
      if ((! subtoolMode) && (! stdinMode) &&
          (! binseqentryMode) && (! smartnetMode)) {
        SeparatorItem (m);
        FormCommandItem (m, "Print", bfp, VIB_MSG_PRINT);
        SeparatorItem (m);
        FormCommandItem (m, "Quit/Q", bfp, VIB_MSG_QUIT);
      } else {
        SeparatorItem (m);
        i = CommandItem (m, "Propagate Top Descriptors", ForcePropagate);
        SetObjectExtra (i, bfp, NULL);
        SeparatorItem (m);
        FormCommandItem (m, "Print", bfp, VIB_MSG_PRINT);
        if (smartnetMode) {
          SeparatorItem (m);
          FormCommandItem (m, "Quit/Q", bfp, VIB_MSG_QUIT);
        }
      }
    } else {
/*#else*/
      if (subtoolMode || stdinMode || binseqentryMode) {
        FormCommandItem (m, "Print", bfp, VIB_MSG_PRINT);
        SeparatorItem (m);
        FormCommandItem (m, "Abort Session", bfp, VIB_MSG_QUIT);
        SeparatorItem (m);
        CommandItem (m, "Accept Changes", SubtoolDoneProc);
      } else if (smartnetMode) {
#ifdef USE_SMARTNET
        FormCommandItem (m, "Print", bfp, VIB_MSG_PRINT);
        SeparatorItem (m);
        FormCommandItem (m, "Abort Session", bfp, VIB_MSG_CLOSE);
        SeparatorItem (m);
        i = CommandItem (m, "Accept Changes", SmartnetDoneProc);
        SetObjectExtra (i, bfp, NULL);
        SeparatorItem (m);
        FormCommandItem (m, "Quit/Q", bfp, VIB_MSG_QUIT);
#endif
      } else {
        SeparatorItem (m);
        FormCommandItem (m, "Print", bfp, VIB_MSG_PRINT);
        SeparatorItem (m);
        FormCommandItem (m, "Quit/Q", bfp, VIB_MSG_QUIT);
      }
    }
/*#endif*/

    m = PulldownMenu (w, "Edit/ E");
    if (subtoolMode || smartnetMode || backupMode) {
      FormCommandItem (m, UNDO_MENU_ITEM, bfp, VIB_MSG_UNDO);
      SeparatorItem (m);
    }
    FormCommandItem (m, COPY_MENU_ITEM, bfp, VIB_MSG_COPY);
    FormCommandItem (m, CLEAR_MENU_ITEM, bfp, VIB_MSG_DELETE);
    SeparatorItem (m);
    if (extraServices) {
      mssgdup = RegisterFormMenuItemName ("SequinDuplicateItem");
      FormCommandItem (m, "Duplicate...", bfp, mssgdup);
      SeparatorItem (m);
    }
    if (genomeCenter != NULL || indexerVersion) {
      SetupEditSecondary (m, bfp);
      SeparatorItem (m);
    }
    mssgseq = RegisterFormMenuItemName ("SequinEditSequenceItem");
    mssgalign = RegisterFormMenuItemName ("SequinEditAlignmentItem");
    mssgsub = RegisterFormMenuItemName ("SequinEditSubmitterItem");
    mssgupd = RegisterFormMenuItemName ("SequinUpdateSeqSubmenu");
    if (indexerVersion)
    {
      mssgupd_idx = RegisterFormMenuItemName ("SequinUpdateSeqSubmenuIndexer");
    }
    mssgext = RegisterFormMenuItemName ("SequinExtendSeqSubmenu");
    mssgfeatprop = RegisterFormMenuItemName ("SequinFeaturePropagate");
    mssgadd = RegisterFormMenuItemName ("SequinAddSeqSubmenu");
    FormCommandItem (m, "Edit Sequence...", bfp, mssgseq);
    FormCommandItem (m, "Alignment Assistant...", bfp, mssgalign);
    FormCommandItem (m, "Edit Submitter Info...", bfp, mssgsub);
    if (indexerVersion) {
      SeparatorItem (m);
      i = CommandItem (m, "Edit History....", EditSequenceHistory);
      SetObjectExtra (i, bfp, NULL);
    }
    SeparatorItem (m);
    
    if (indexerVersion)
    {
      /* indexer version */
      sub = SubMenu (m, "Indexer Update Sequence");
      SetFormMenuItem (bfp, mssgupd_idx, (IteM) sub);
      i = CommandItem (sub, "Single Sequence", TestUpdateSequenceIndexer);
      SetObjectExtra (i, bfp, NULL);
      if (useEntrez) {
        i = CommandItem (sub, "Download Accession", UpdateSequenceViaDownloadIndexer);
        SetObjectExtra (i, bfp, NULL);    
      }
      i = CommandItem (sub, "Multiple Sequences", TestUpdateSequenceSetIndexer);
      SetObjectExtra (i, bfp, NULL);
      
      /* public version */
      sub = SubMenu (m, "Public Update Sequence");
      i = CommandItem (sub, "Single Sequence", TestUpdateSequenceSubmitter);
      SetObjectExtra (i, bfp, NULL);
      if (useEntrez) {
        i = CommandItem (sub, "Download Accession", UpdateSequenceViaDownloadSubmitter);
        SetObjectExtra (i, bfp, NULL);    
      }
      i = CommandItem (sub, "Multiple Sequences", TestUpdateSequenceSetSubmitter);
      SetObjectExtra (i, bfp, NULL);
      
      SeparatorItem (m);
    }
    else
    {
      sub = SubMenu (m, "Update Sequence");
      SetFormMenuItem (bfp, mssgupd, (IteM) sub);
      i = CommandItem (sub, "Single Sequence", TestUpdateSequenceSubmitter);
      SetObjectExtra (i, bfp, NULL);
      if (useEntrez) {
        i = CommandItem (sub, "Download Accession", UpdateSequenceViaDownloadSubmitter);
        SetObjectExtra (i, bfp, NULL);    
      }
      i = CommandItem (sub, "Multiple Sequences", TestUpdateSequenceSetSubmitter);
      SetObjectExtra (i, bfp, NULL);
    }
    
    sub = SubMenu (m, "Extend Sequence");
    SetFormMenuItem (bfp, mssgext, (IteM) sub);
    i = CommandItem (sub, "Read FASTA File...", ExtendSeqWithFASTA);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (sub, "Read Sequence Record...", ExtendSeqWithRec);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (sub, "Read FASTA or ASN.1 Set", ExtendFastaSet);
    SetObjectExtra (i, bfp, NULL);
    if (useEntrez) {
      i = CommandItem (sub, "Download Accession...", ExtendSeqWithAcc);
      SetObjectExtra (i, bfp, NULL);
    }
            
    SeparatorItem (m);
    i = CommandItem (m, "Feature Propagate...", NewFeaturePropagate);
    SetObjectExtra (i, bfp, NULL);
    SetFormMenuItem (bfp, mssgfeatprop, i);
    SeparatorItem (m);
    sub = SubMenu (m, "Add Sequence");
    SetFormMenuItem (bfp, mssgadd, (IteM) sub);
    i = CommandItem (sub, "Add FASTA File...", AddSeqWithFASTA);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (sub, "Add Sequence Record...", AddSeqWithRec);
    SetObjectExtra (i, bfp, NULL);
    if (! extraServices) {
      SeparatorItem (m);
      i = CommandItem (m, "Parse File to Source", ParseFileToSource);
      SetObjectExtra (i, bfp, NULL);
    }

    m = PulldownMenu (w, "Search/ R");
    i = CommandItem (m, "Find ASN.1.../ F", FindStringProc);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (m, "Find FlatFile.../ G", FindFlatfileProc);
    SetObjectExtra (i, bfp, NULL);
    SeparatorItem (m);
    i = CommandItem (m, "Find by Gene...", FindGeneProc);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (m, "Find by Protein...", FindProtProc);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (m, "Find by Position...", FindPosProc);
    SetObjectExtra (i, bfp, NULL);
    SeparatorItem (m);
    if (indexerVersion) {
      sub = SubMenu (m, "Validate/ V");
      i = CommandItem (sub, "Validate Record/ R", ValSeqEntryProc);
      SetObjectExtra (i, bfp, NULL);
      i = CommandItem (sub, "Validate no Alignments/ A", ValSeqEntryProcNoAln);
      SetObjectExtra (i, bfp, NULL);
      i = CommandItem (sub, "Validate check Inference", ValSeqEntryProcInfAccn);
      SetObjectExtra (i, bfp, NULL);
      SeparatorItem (sub);
      i = CommandItem (sub, "Validate Inst", ValSeqEntryProcInst);
      SetObjectExtra (i, bfp, NULL);
      i = CommandItem (sub, "Validate Hist", ValSeqEntryProcHist);
      SetObjectExtra (i, bfp, NULL);
      i = CommandItem (sub, "Validate Context", ValSeqEntryProcContext);
      SetObjectExtra (i, bfp, NULL);
      i = CommandItem (sub, "Validate Graph", ValSeqEntryProcGraph);
      SetObjectExtra (i, bfp, NULL);
      i = CommandItem (sub, "Validate Set", ValSeqEntryProcSet);
      SetObjectExtra (i, bfp, NULL);
      i = CommandItem (sub, "Validate Feat", ValSeqEntryProcFeat);
      SetObjectExtra (i, bfp, NULL);
      i = CommandItem (sub, "Validate Desc", ValSeqEntryProcDesc);
      SetObjectExtra (i, bfp, NULL);
    } else {
      i = CommandItem (m, "Validate/ V", ValSeqEntryProc);
      SetObjectExtra (i, bfp, NULL);
    }
#ifdef USE_SPELL
    SeparatorItem (m);
    i = CommandItem (m, "Spell Check...", SpellCheckSeqEntryProc);
    SetObjectExtra (i, bfp, NULL);
#endif
/*#ifdef USE_BLAST*/
    if (useBlast) {
      SeparatorItem (m);
      if (indexerVersion) {
        i = CommandItem (m, "CDD BLAST...", SimpleCDDBlastProc);
        SetObjectExtra (i, bfp, NULL);
      }
      if (extraServices) {
        sub = SubMenu (m, "CDD Search");
        i = CommandItem (sub, "Features", SimpleCDDSearchFeatProc);
        SetObjectExtra (i, bfp, NULL);
        i = CommandItem (sub, "Alignments", SimpleCDDSearchAlignProc);
        SetObjectExtra (i, bfp, NULL);
      } else {
        i = CommandItem (m, "CDD Search", SimpleCDDSearchFeatProc);
        SetObjectExtra (i, bfp, NULL);
      }
    }
/*#endif*/
    SeparatorItem (m);
    sub = SubMenu (m, "Vector Screen");
    i = CommandItem (sub, "UniVec", SimpleUniVecScreenProc);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (sub, "UniVec Core", SimpleUniVecCoreScreenProc);
    SetObjectExtra (i, bfp, NULL);
    SeparatorItem (m);
    i = CommandItem (m, "ORF Finder...", FindOrf);
    SetObjectExtra (i, bfp, NULL);
    /*
    i = CommandItem (m, "Repeat Finder...", FindAlu);
    SetObjectExtra (i, bfp, NULL);
    */
    SeparatorItem (m);
    i = CommandItem (m, "Select Target...", DoChangeTarget);
    SetObjectExtra (i, bfp, NULL);

    /*
    if (! indexerVersion) {
      m = PulldownMenu (w, "Options");
      sub = SubMenu (m, "Font Selection");
      i = CommandItem (sub, "Display Font...", DisplayFontChangeProc);
      SetObjectExtra (i, bfp, NULL);
      SeparatorItem (m);
      CreateLegendItem (m, bfp);
    }
    */

/*#ifdef EXTRA_SERVICES*/
    if (extraServices) {
      m = PulldownMenu (w, "Special/ S");
      SetupSpecialMenu (m, bfp);
      m = PulldownMenu (w, "Projects");
      MakeSpecialProjectsMenu (m, bfp);
    }
/*#endif*/

    m = PulldownMenu (w, "Misc");
    /*
    CommandItem (m, "Style Manager...", StyleManagerProc);
    SeparatorItem (m);
    */
    CommandItem (m, "Net Configure...", NetConfigureProc);
    if (useEntrez) {
      /*
      SeparatorItem (m);
      CommandItem (m, "Entrez2 Query...", Entrez2QueryProc);
      */
/*
#ifndef WIN16
      if (BiostrucAvail ()) {
        SeparatorItem (m);
        CommandItem (m, "Cn3D Window...", Cn3DWinShowProc);
      }
#endif
*/
    }
    if (useDesktop) {
      SeparatorItem (m);
      VSMAddToMenu (m, VSM_DESKTOP);
    }

    CreateAnalysisMenu (w, bfp, TRUE, FALSE);

    m = PulldownMenu (w, "Annotate/ A");
    SetupNewFeaturesMenu (m, bfp);
    SeparatorItem (m);
    sub = SubMenu (m, "Batch Feature Apply");
    SetupBatchApplyMenu (sub, bfp);
    sub = SubMenu (m, "Batch Feature Edit");
    SetupBatchEditMenu (sub, bfp);
    i = CommandItem (m, "Batch Apply Molecule Type", ExternalApplyMoleculeType);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (m, "Set Release Date", SetReleaseDate);
    SetObjectExtra (i, bfp, NULL);
    SeparatorItem (m);
    i = CommandItem (m, "ORF Finder", FindOrf);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (m, "Import Source Table", ParseFileToSource);
    SetObjectExtra (i, bfp, NULL);
    SeparatorItem (m);
    sub = SubMenu (m, "Publications");
    SetupNewPublicationsMenu (sub, bfp);
    SeparatorItem (m);
    sub = SubMenu (m, "Descriptors");
    SetupNewDescriptorsMenu (sub, bfp);
    SeparatorItem (m);
    sub = SubMenu (m, "Advanced Table Readers");
    i = CommandItem (sub, "Load Structured Comments from Table", SubmitterCreateStructuredComments);
    SetObjectExtra (i, bfp, NULL);
    i = CommandItem (m, "Sort Unique Count By Group", SUCSubmitterProc);
    SetObjectExtra (i, bfp, NULL);

    if (indexerVersion) {
      m = PulldownMenu (w, "Options");
      sub = SubMenu (m, "Font Selection");
      i = CommandItem (sub, "Display Font...", DisplayFontChangeProc);
      SetObjectExtra (i, bfp, NULL);
      /*
      SeparatorItem (m);
      CreateLegendItem (m, bfp);
      */
      SeparatorItem (m);
      sub = SubMenu (m, "Layout Override");
      CreateNewLayoutMenu (sub, bfp);
    }
  }
}
#endif

static void GetRidCallback (SeqEntryPtr sep, Pointer mydata, Int4 index, Int2 indent)

{
  BioseqPtr     bsp;
  BioseqSetPtr  bssp;
  SeqAnnotPtr   nextsap;
  Pointer PNTR  prevsap;
  SeqAnnotPtr   sap;

  if (sep == NULL || sep->data.ptrvalue == NULL) return;
  if (IS_Bioseq (sep)) {
    bsp = (BioseqPtr) sep->data.ptrvalue;
    sap = bsp->annot;
    prevsap = (Pointer PNTR) &(bsp->annot);
  } else if (IS_Bioseq_set (sep)) {
    bssp = (BioseqSetPtr) sep->data.ptrvalue;
    sap = bssp->annot;
    prevsap = (Pointer PNTR) &(bssp->annot);
  } else return;
  while (sap != NULL) {
    nextsap = sap->next;
    if (sap->data == NULL) {
      *(prevsap) = sap->next;
      sap->next = NULL;
      SeqAnnotFree (sap);
    } else {
      prevsap = (Pointer PNTR) &(sap->next);
    }
    sap = nextsap;
  }
}

static void GetRidOfEmptyAnnotTables (Uint2 entityID)

{
  SeqEntryPtr  sep;

  if (entityID < 1) return;
  sep = GetTopSeqEntryForEntityID (entityID);
  if (sep == NULL) return;
  SeqEntryExplore (sep, NULL, GetRidCallback);
}

static CharPtr deleteProtMsg =
"The protein product of a CDS (shown in the /translation qualifier)\n\
is actually a separate data element in the record.  Unless explicitly\n\
deleted, it will remain hidden in the record after you delete the CDS.";

static CharPtr deleteGeneMsg =
"The /gene qualifier is generated from an overlapping gene feature.\n\
If you delete a CDS you may also want to delete this separate gene.";

static CharPtr deleteCdnaMsg =
"The cDNA product of an mRNA is actually a separate data element\n\
in the record.  Unless explicitly deleted, it will remain hidden\n\
in the record after you delete the mRNA.";

typedef struct deletecdsoptions {
  Boolean delete_feature;
  Boolean delete_gene;
  Boolean delete_cdna;
  Boolean delete_protein;
} DeleteCDSOptionsData, PNTR DeleteCDSOptionsPtr;

static void GetDeleteCDSOptions (DeleteCDSOptionsPtr dcop)
{
  WindoW w;
  GrouP  h, g, c;
  ButtoN b;
  ModalAcceptCancelData acd;
  ButtoN delete_gene = NULL;
  ButtoN delete_cdna = NULL;
  ButtoN delete_protein = NULL;
  
  if (dcop == NULL) return;
  if (!dcop->delete_gene
      && !dcop->delete_cdna
      && !dcop->delete_protein) {
    return;
  }
  
  acd.accepted = FALSE;
  acd.cancelled = FALSE;
  
  w = ModalWindow(-20, -13, -10, -10, NULL);
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  g = HiddenGroup (h, 2, 0, NULL);
  SetGroupSpacing (g, 10, 10);
  if (dcop->delete_gene) {
    delete_gene = CheckBox (g, "Delete overlapping gene", NULL);
    MultiLinePrompt (g, deleteGeneMsg, 30 * stdCharWidth, systemFont);
  }
  
  if (dcop->delete_cdna) {
    delete_cdna = CheckBox (g, "Delete cDNA", NULL);
    SetStatus (delete_cdna, TRUE);
    MultiLinePrompt (g, deleteCdnaMsg, 30 * stdCharWidth, systemFont);
  }

  if (dcop->delete_protein) {
    delete_protein = CheckBox (g, "Delete protein product", NULL);
    SetStatus (delete_protein, TRUE);
    MultiLinePrompt (g, deleteProtMsg, 30 * stdCharWidth, systemFont);
  }

  c = HiddenGroup (h, 3, 0, NULL);
  SetGroupSpacing (c, 10, 10);
  b = PushButton (c, "Accept", ModalAcceptButton);
  SetObjectExtra (b, &acd, NULL);
  b = PushButton (c, "Cancel", ModalCancelButton);
  SetObjectExtra (b, &acd, NULL);
  AlignObjects (ALIGN_CENTER, (HANDLE) g, (HANDLE) c, NULL);
  
  Show(w); 
  Select (w);
  while (!acd.accepted && ! acd.cancelled)
  {
    ProcessExternalEvent ();
    Update ();
  }
  ProcessAnEvent ();
  Remove (w);
  if (acd.accepted)
  {
    dcop->delete_feature = TRUE;
    dcop->delete_gene = delete_gene == NULL ? FALSE : GetStatus (delete_gene);
    dcop->delete_cdna = delete_cdna == NULL ? FALSE : GetStatus (delete_cdna);
    dcop->delete_protein = delete_protein == NULL ? FALSE : GetStatus (delete_protein);
  }
  else
  {
    dcop->delete_feature = FALSE;
    dcop->delete_gene = FALSE;
    dcop->delete_cdna = FALSE;
    dcop->delete_protein = FALSE;
  }
}

static Boolean DeleteSelectedFeatureOrDescriptor (GatherContextPtr gcp)

{
  BioseqPtr       bsp;
  BioseqPtr       cdna;
  Uint2           entityID;
  SeqFeatPtr      gene;
  Uint4           itemID;
  Uint2           itemtype;
  BioseqSetPtr    nps;
  BioseqSetPtr    parent;
  SeqFeatPtr      sfp;
  SeqIdPtr        sip;
  DeleteCDSOptionsData dcod;
  SeqMgrFeatContext    fcontext;
  ObjValNodePtr        ovn;
  SeqDescrPtr          sdp;
#ifdef USE_SMARTNET
  UserObjectPtr   uop;
#endif

  if (gcp->thistype != OBJ_SEQDESC && gcp->thistype != OBJ_SEQFEAT) {
    /* This function should only handle descriptors and features */
    return FALSE;
  } else if (gcp->thisitem == NULL) {
    return TRUE;
  }

#ifdef USE_SMARTNET
  /* This code will prevent from deletion SMART User Object */

  if(gcp->thistype == OBJ_SEQDESC) {
      if((sdp = (SeqDescrPtr) gcp->thisitem) != NULL) {
          if(sdp->choice == 14 &&
             ((uop = ( UserObjectPtr) sdp->data.ptrvalue) != NULL)) {
              if(!StringCmp(uop->_class, SMART_OBJECT_LABEL)) {
                  Message(MSG_ERROR, "You may not delete SMART Object Label");
                  return TRUE;
              }
          }
      }
  }
#endif

  /* delete the descriptor */
  if (gcp->thistype == OBJ_SEQDESC) {
    sdp = (SeqDescrPtr) gcp->thisitem;
    if (sdp != NULL && sdp->extended != 0) {
      ovn = (ObjValNodePtr) sdp;
      ovn->idx.deleteme = TRUE;
      return TRUE;
    } else {
      return FALSE;
    }
  }

  /* because it wasn't a descriptor, this must be a feature */
  sfp = (SeqFeatPtr) gcp->thisitem;
  bsp = NULL;
  gene = NULL;
  cdna = NULL;
  nps = NULL;
  entityID = gcp->entityID;
  itemID = gcp->itemID;
  itemtype = gcp->thistype;

  /* When deleting some features, there are other objects that the user may wish to
   * remove at the same time.
   * When removing a coding region, the user may wish to also do the following:
   *   1) Remove the overlapping gene
   *   2) Remove the product protein
   * When removing an mRNA, the user may wish to also do the following:
   *   1) Remove the cDNA product
   *   2) Remove the nucprotset parent
   */
  if (sfp->idx.subtype == FEATDEF_mRNA && sfp->product != NULL) {
    sip = SeqLocId (sfp->product);
    if (sip != NULL) {
      cdna = BioseqFind (sip);
      if (cdna != NULL) {
        if (cdna->idx.parenttype == OBJ_BIOSEQSET) {
          parent = (BioseqSetPtr) cdna->idx.parentptr;
          while (parent != NULL) {
            if (parent->_class == BioseqseqSet_class_nuc_prot) {
              nps = parent;
            }
            if (parent->idx.parenttype == OBJ_BIOSEQSET) {
              parent = (BioseqSetPtr) parent->idx.parentptr;
            } else {
              parent = NULL;
            }
          }
        }
      }
    }
  } else if (sfp->idx.subtype == FEATDEF_CDS) {
    if (sfp->product != NULL) {
      sip = SeqLocId (sfp->product);
      if (sip != NULL) {
        bsp = BioseqFind (sip);
      }
    }
    if (SeqMgrGetGeneXref(sfp) == NULL) {
      gene = SeqMgrGetOverlappingGene (sfp->location, &fcontext);
    }
  }

  MemSet (&dcod, 0, sizeof (dcod));
  if (sfp != NULL) {
    dcod.delete_feature = TRUE;
  }
  if (nps != NULL || cdna != NULL) {
    dcod.delete_cdna = TRUE;
  }
  if (bsp != NULL) {
    dcod.delete_protein = TRUE;
  }

  if (gene != NULL) {
    dcod.delete_gene = TRUE;
  }

  GetDeleteCDSOptions (&dcod);
  if (!dcod.delete_feature) {
    return TRUE;
  }


  /* delete the feature */
  sfp->idx.deleteme = TRUE;

  /* delete cdna product */
  if (dcod.delete_cdna) {
    if (nps != NULL) {
      nps->idx.deleteme = TRUE;
    } else if (cdna != NULL) {
      cdna->idx.deleteme = TRUE;
    }
  }

  /* delete protein product */
  if (bsp != NULL && dcod.delete_protein) {
    bsp->idx.deleteme = TRUE;
  }

  /* delete overlapping gene */
  if (gene != NULL && dcod.delete_gene) {
    gene->idx.deleteme = TRUE;
  }

  /* Note - delete marked objects is called by the calling function */

  return TRUE;
}

static void TermSelectionFormMessage (ForM f, Int2 mssg)

{
  BaseFormPtr  bfp;

  bfp = (BaseFormPtr) GetObjectExtra (f);
  if (bfp != NULL) {
    switch (mssg) {
      case VIB_MSG_CLOSE :
        Hide (f);
        break;
      case VIB_MSG_QUIT :
        QuitProc ();
        break;
      default :
        break;
    }
  }
}

static void DocumentSummaryFormMessage (ForM f, Int2 mssg)

{
  BaseFormPtr  bfp;

  bfp = (BaseFormPtr) GetObjectExtra (f);
  if (bfp != NULL) {
    switch (mssg) {
      case VIB_MSG_CHANGE :
        EnableAnalysisItems (bfp, TRUE);
        break;
      case VIB_MSG_CLOSE :
        Hide (f);
        break;
      case VIB_MSG_QUIT :
        QuitProc ();
        break;
      default :
        break;
    }
  }
}

static void SequinMedlineFormMessage (ForM f, Int2 mssg)

{
  BaseFormPtr  bfp;

  bfp = (BaseFormPtr) GetObjectExtra (f);
  if (bfp != NULL) {
    switch (mssg) {
      case VIB_MSG_CLOSE :
        Remove (f);
        break;
      case VIB_MSG_QUIT :
        QuitProc ();
        break;
      default :
        break;
    }
  }
}

typedef struct seltbl {
  size_t             count;
  SelStructPtr PNTR  selarray;
  Boolean            geneasked;
  Boolean            productasked;
  Boolean            removegene;
  Boolean            removeproduct;
} SelTbl, PNTR SelTblPtr;

static Boolean MarkSelectedItem (GatherObjectPtr gop)

{
  Int2           L, R, mid;
  SelStructPtr   ssp;
  SelTblPtr      tbl;
  MsgAnswer      ans;
  BioseqPtr      bsp;
  ObjValNodePtr  ovn;
  SeqAlignPtr    sap;
  SeqDescrPtr    sdp;
  SeqFeatPtr     sfp, gene;
  SeqGraphPtr    sgp;

  if (gop == NULL) return TRUE;
  tbl = (SelTblPtr) gop->userdata;
  if (tbl == NULL) return TRUE;

  L = 0;
  R = tbl->count - 1;
  while (L <= R) {
    mid = (L + R) / 2;
    ssp = tbl->selarray [mid];
    if (ssp == NULL) return TRUE;
    if (ssp->entityID > gop->entityID) {
      R = mid - 1;
    } else if (ssp->entityID < gop->entityID) {
      L = mid + 1;
    } else if (ssp->itemtype > gop->itemtype) {
      R = mid - 1;
    } else if (ssp->itemtype < gop->itemtype) {
      L = mid + 1;
    } else if (ssp->itemID > gop->itemID) {
      R = mid - 1;
    } else if (ssp->itemID < gop->itemID) {
      L = mid + 1;
    } else if (gop->dataptr != NULL &&
               ssp->entityID == gop->entityID &&
               ssp->itemtype == gop->itemtype &&
               ssp->itemID == gop->itemID) {
      switch (gop->itemtype) {
        case OBJ_SEQDESC :
          sdp = (SeqDescrPtr) gop->dataptr;
          if (sdp != NULL && sdp->extended != 0) {
            ovn = (ObjValNodePtr) sdp;
            ovn->idx.deleteme = TRUE;
          }
          break;
        case OBJ_SEQFEAT :
          sfp = (SeqFeatPtr) gop->dataptr;
          if (sfp != NULL) {
            sfp->idx.deleteme = TRUE;
            if (SeqMgrGetGeneXref (sfp) == NULL) {
              gene = SeqMgrGetOverlappingGene (sfp->location, NULL);
              if (gene != NULL) {
                if (! tbl->geneasked) {
                  ans = Message (MSG_YN, "Remove overlapping gene?");
                  if (ans == ANS_YES) {
                    tbl->removegene = TRUE;
                  }
                  tbl->geneasked = TRUE;
                }
                if (tbl->removegene) {
                  gene->idx.deleteme = TRUE;
                }
              }
            }
            if (sfp->data.choice == SEQFEAT_CDREGION) {
              bsp = BioseqFind (SeqLocId (sfp->product));
              if (bsp != NULL) {
                if (! tbl->productasked) {
                  ans = Message (MSG_YN, "Remove protein products?");
                  if (ans == ANS_YES) {
                    tbl->removeproduct = TRUE;
                  }
                  tbl->productasked = TRUE;
                }
                if (tbl->removeproduct) {
                  bsp->idx.deleteme = TRUE;
                }
              }
            }
          }
          break;
        case OBJ_SEQALIGN :
          sap = (SeqAlignPtr) gop->dataptr;
          if (sap != NULL) {
            sap->idx.deleteme = TRUE;
          }
          break;
        case OBJ_SEQGRAPH :
          sgp = (SeqGraphPtr) gop->dataptr;
          if (sgp != NULL) {
            sgp->idx.deleteme = TRUE;
          }
          break;
        default :
          break;
      }
      return TRUE;
    }
  }

  return TRUE;
}

static int LIBCALLBACK SortSelStructByIDs (VoidPtr ptr1, VoidPtr ptr2)

{
  SelStructPtr  ssp1, ssp2;

  if (ptr1 == NULL || ptr2 == NULL) return 0;
  ssp1 = *((SelStructPtr PNTR) ptr1);
  ssp2 = *((SelStructPtr PNTR) ptr2);
  if (ssp1 == NULL || ssp2 == NULL) return 0;

  if (ssp1->entityID > ssp2->entityID) return 1;
  if (ssp1->entityID < ssp2->entityID) return -1;

  if (ssp1->itemtype > ssp2->itemtype) return 1;
  if (ssp1->itemtype < ssp2->itemtype) return -1;

  if (ssp1->itemID > ssp2->itemID) return 1;
  if (ssp1->itemID < ssp2->itemID) return -1;

  return 0;
}

static void DeleteMultipleSelections (Uint2 entityID, SelStructPtr selhead)

{
  size_t             count;
  SelStructPtr       sel;
  SelStructPtr PNTR  selarray;
  SelTbl             tbl;

  if (entityID < 1 || selhead == NULL) return;

  for (count = 0, sel = selhead; sel != NULL; sel = sel->next, count++) continue;
  if (count < 1) return;

  selarray = (SelStructPtr PNTR) MemNew (sizeof (SelStructPtr) * (size_t) (count + 1));
  if (selarray == NULL) return;
  for (count = 0, sel = selhead; sel != NULL; sel = sel->next, count++) {
    selarray [count] = sel;
  }

  HeapSort (selarray, (size_t) count, sizeof (SelStructPtr), SortSelStructByIDs);

  MemSet ((Pointer) &tbl, 0, sizeof (tbl));
  tbl.count = count;
  tbl.selarray = selarray;

  GatherObjectsInEntity (entityID, 0, NULL, MarkSelectedItem, (Pointer) &tbl, NULL);

  MemFree (selarray);
  DeleteMarkedObjects (entityID, 0, NULL);
}

extern void SequinSeqViewFormMessage (ForM f, Int2 mssg)

{
  MsgAnswer     ans;
  BaseFormPtr   bfp;
  BioseqPtr     bsp;
  Uint4         entityID;
  Uint4         itemID;
  Int2          mssgalign;
  Int2          mssgdup;
  Int2          mssgseq;
  Int2          mssgsub;
  SeqEntryPtr   oldsep;
  OMProcControl  ompc;
/******** COLOMBE
  SeqAlignPtr   salp; ********COLOMBE END*/
  SeqAnnotPtr   sap;
  SelStructPtr  sel;
  SeqEntryPtr   sep;
  SeqIdPtr      sip;

  bfp = (BaseFormPtr) GetObjectExtra (f);
  if (bfp != NULL) {
    switch (mssg) {
      case VIB_MSG_SAVE :
        if (FixSpecialCharacters (bfp->input_entityID))
        {
          SaveSeqSubmitProc (bfp, FALSE);
        }
        break;
      case VIB_MSG_SAVE_AS :
        if (FixSpecialCharacters (bfp->input_entityID))
        {
          SaveSeqSubmitProc (bfp, TRUE);
        }
        break;
      case VIB_MSG_CLOSE :
        CloseProc (bfp);
        break;
      case VIB_MSG_QUIT :
        QuitProc ();
        break;
      case VIB_MSG_ACCEPT :
        ProcessDoneButton (f);
        break;
      case VIB_MSG_RESET :
#ifdef USE_SMARTNET
        SmartResetProc ((IteM)f);
#endif
        break;
      case VIB_MSG_IMPORT :
        /* ReadFastaProc (); */
        break;
      case VIB_MSG_CUT :
        break;
      case VIB_MSG_COPY :
        break;
      case VIB_MSG_UNDO :
        if (subtoolMode || smartnetMode || backupMode) {
          if (FileLength (SEQUIN_EDIT_PREV_FILE) > 0) {
            sep = GetTopSeqEntryForEntityID (subtoolEntityID);
            if (Message (MSG_YN, "Restore from backup?") == ANS_YES) {
              oldsep = RestoreFromFile (SEQUIN_EDIT_PREV_FILE);
              ReplaceSeqEntryWithSeqEntry (sep, oldsep, TRUE);
              subtoolEntityID = ObjMgrGetEntityIDForChoice (sep);
              ObjMgrSetDirtyFlag (subtoolEntityID, TRUE);
              ObjMgrSendMsg (OM_MSG_UPDATE, subtoolEntityID, 0, 0);
            }
          }
        }
        break;
      case VIB_MSG_PASTE :
        break;
      case VIB_MSG_DELETE :
        sel = ObjMgrGetSelected ();
        if (sel != NULL) {
          if (sel->itemtype == OBJ_BIOSEQ && sel->next == NULL) {
            if (! indexerVersion) {
              ans = Message (MSG_OKC, "Are you sure you want to delete this Bioseq?");
              if (ans == ANS_CANCEL) return;
              ans = Message (MSG_OKC, "Are you REALLY sure you want to delete this Bioseq?");
              if (ans == ANS_CANCEL) return;
            }
            bsp = GetBioseqGivenIDs (sel->entityID, sel->itemID, sel->itemtype);
            sip = SeqIdDup (bsp->id);
            entityID = sel->entityID;
            MemSet ((Pointer) (&ompc), 0, sizeof (OMProcControl));
            ompc.do_not_reload_from_cache = TRUE;
            ompc.input_entityID = sel->entityID;
            ompc.input_itemID = sel->itemID;
            ompc.input_itemtype = sel->itemtype;
            if (! DetachDataForProc (&ompc, FALSE)) {
              Message (MSG_ERROR, "DetachDataForProc failed");
            }
            sep = GetTopSeqEntryForEntityID (entityID);
            SeqAlignBioseqDeleteByIdFromSeqEntry (sep, sip);
            SeqIdFree (sip);
            /* see VSeqMgrDeleteProc - probably leaving one dangling SeqEntry */
            BioseqFree (bsp);
            ObjMgrSetDirtyFlag (entityID, TRUE);
            ObjMgrSendMsg (OM_MSG_UPDATE, entityID, 0, 0);
            ObjMgrDeSelect (0, 0, 0, 0, NULL);
            Update ();
          }
          if (sel->itemtype == OBJ_SEQANNOT && sel->next == NULL && extraServices) {
            entityID = sel->entityID;
            MemSet ((Pointer) (&ompc), 0, sizeof (OMProcControl));
            ompc.do_not_reload_from_cache = TRUE;
            ompc.input_entityID = sel->entityID;
            ompc.input_itemID = sel->itemID;
            ompc.input_itemtype = sel->itemtype;
            if (! DetachDataForProc (&ompc, FALSE)) {
              Message (MSG_ERROR, "DetachDataForProc failed");
            }
            ObjMgrSetDirtyFlag (entityID, TRUE);
            ObjMgrSendMsg (OM_MSG_UPDATE, entityID, 0, 0);
            ObjMgrDeSelect (0, 0, 0, 0, NULL);
            Update ();
          }
          if (sel->next == NULL && (sel->itemtype == OBJ_SEQDESC || sel->itemtype == OBJ_SEQFEAT)) {
            entityID = sel->entityID;
            GatherItem (sel->entityID, sel->itemID, sel->itemtype,
                        NULL, DeleteSelectedFeatureOrDescriptor);
            DeleteMarkedObjects (entityID, 0, NULL);
            GetRidOfEmptyAnnotTables (entityID);
            sep = GetTopSeqEntryForEntityID (entityID);
            RenormalizeNucProtSets (sep, TRUE);
            ObjMgrSetDirtyFlag (entityID, TRUE);
            ObjMgrSendMsg (OM_MSG_UPDATE, entityID, 0, 0);
            ObjMgrDeSelect (0, 0, 0, 0, NULL);
            Update ();
          } else {
            /* Message (MSG_OK, "Unable to delete multiple objects"); */
            ans = Message (MSG_YN, "Are you sure you want to delete multiple objects?");
            if (ans == ANS_NO) return;
            entityID = bfp->input_entityID;
            DeleteMultipleSelections (entityID, sel);
            GetRidOfEmptyAnnotTables (entityID);
            sep = GetTopSeqEntryForEntityID (entityID);
            RenormalizeNucProtSets (sep, TRUE);
            ObjMgrSetDirtyFlag (entityID, TRUE);
            ObjMgrSendMsg (OM_MSG_UPDATE, entityID, 0, 0);
            ObjMgrDeSelect (0, 0, 0, 0, NULL);
          }
        } else {
          Message (MSG_OK, "Nothing selected");
        }
        break;
      case VIB_MSG_CHANGE :
        EnableFeaturesPerTarget (bfp);
        EnableAnalysisItems (bfp, FALSE);
        EnableEditSeqAlignAndSubItems (bfp);
        break;
      case VIB_MSG_SELECT :
        EnableEditAlignItem (bfp);
        break;
      default :
        mssgseq = RegisterFormMenuItemName ("SequinEditSequenceItem");
        mssgalign = RegisterFormMenuItemName ("SequinEditAlignmentItem");
        mssgsub = RegisterFormMenuItemName ("SequinEditSubmitterItem");
        mssgdup = RegisterFormMenuItemName ("SequinDuplicateItem");
        if (mssg == mssgseq) {
          bsp =  GetBioseqGivenIDs (bfp->input_entityID, bfp->input_itemID, bfp->input_itemtype);
          if (bsp != NULL) {
            WatchCursor ();
            Update ();
            SeqEntrySetScope (NULL);
            GatherProcLaunch (OMPROC_EDIT, FALSE, bfp->input_entityID, bfp->input_itemID,
                              bfp->input_itemtype, 0, 0, bfp->input_itemtype, 0);
            ArrowCursor ();
            Update ();
          }
        } else if (mssg == mssgalign) {
          sel = ObjMgrGetSelected ();
          if (sel != NULL &&
              (sel->itemtype == OBJ_SEQALIGN || sel->itemtype == OBJ_SEQHIST_ALIGN)) {
            WatchCursor ();
            Update ();
            SeqEntrySetScope (NULL);
            GatherProcLaunch (OMPROC_EDIT, FALSE, sel->entityID, sel->itemID,
                              sel->itemtype, 0, 0, sel->itemtype, 0);
            ArrowCursor ();
            Update ();
          } else /* if (sel == NULL) */ {
            sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
/******** COLOMBE
            salp = FindSeqAlignInSeqEntry (sep, OBJ_SEQALIGN);
            if (salp != NULL) {
              itemID = GetItemIDGivenPointer (bfp->input_entityID, OBJ_SEQALIGN, (Pointer) salp);
              if (itemID != 0) {
                WatchCursor ();
                Update ();
                GatherProcLaunch (OMPROC_EDIT, FALSE, bfp->input_entityID, itemID,
                                  OBJ_SEQALIGN, 0, 0, OBJ_SEQALIGN, 0);
                ArrowCursor ();
                Update ();
              }
            }
COLOMBE END *******************/
            sap  = (SeqAnnotPtr)FindSeqAlignInSeqEntry (sep, OBJ_SEQANNOT);
            if (sap) {
              itemID = GetItemIDGivenPointer (bfp->input_entityID, OBJ_SEQANNOT, (Pointer)sap);
              if (itemID != 0) {
                WatchCursor ();
                Update ();
                SeqEntrySetScope (NULL);
                GatherProcLaunch (OMPROC_EDIT, FALSE, bfp->input_entityID, itemID, OBJ_SEQANNOT, 0, 0, OBJ_SEQANNOT, 0);
                ArrowCursor ();
                Update ();
              }
            }
          }
        } else if (mssg == mssgsub) {
          EditSubmitBlock (bfp);
        } else if (mssg == mssgdup) {
          sel = ObjMgrGetSelected ();
          if (sel != NULL &&
              (sel->itemtype == OBJ_SEQDESC || sel->itemtype == OBJ_SEQFEAT)) {
            stdedprocs.duplicateExisting = TRUE;
            WatchCursor ();
            Update ();
            SeqEntrySetScope (NULL);
            GatherProcLaunch (OMPROC_EDIT, FALSE, sel->entityID, sel->itemID,
                              sel->itemtype, 0, 0, sel->itemtype, 0);
            ArrowCursor ();
            Update ();
            stdedprocs.duplicateExisting = FALSE;
          }
        }
        break;
    }
  }
}

static void SequinSeqEditFormMessage (ForM f, Int2 mssg)

{
  BaseFormPtr  bfp;

  bfp = (BaseFormPtr) GetObjectExtra (f);
  if (bfp != NULL) {
    switch (mssg) {
      case VIB_MSG_QUIT :
        QuitProc ();
        break;
      default :
        break;
    }
  }
}

/* from salfiles.c */
static Int4 AccessionToGi (CharPtr string, CharPtr program)
{
   /*
   CharPtr str;
   LinkSetPtr lsp;
   Int4 gi;

   EntrezInit (program, TRUE, NULL);
   str = MemNew (StringLen (string) + 10);
   sprintf (str, "\"%s\" [ACCN]", string);
   lsp = EntrezTLEvalString (str, TYP_NT, -1, NULL, NULL);
   MemFree (str);
   if (lsp == NULL)
       return 0;
   if (lsp->num <= 0) {
       LinkSetFree (lsp);
       return 0;
   }
   gi = lsp->uids [0];
   LinkSetFree (lsp);
   EntrezFini ();
   return gi;
   */
   return 0;
}

static SeqEntryPtr LIBCALLBACK SeqEdDownload (CharPtr program, CharPtr accession,
                                              Int4 uid, Boolean is_na, BoolPtr is_new)

{
  BioseqPtr    bsp;
  SeqEntryPtr  sep = NULL;
  SeqId        sid;
  /*
  LinkSetPtr   lsp;
  Int2         seqtype;
  Char         str [64];
  */

  if (is_new != NULL) {
     *is_new = TRUE;
  }
/*********/
  if (uid==0 && ! StringHasNoText (accession))
    uid = AccessionToGi (accession, program);
  if (uid > 0) {
    sid.choice = SEQID_GI;
    sid.data.intvalue = uid;
    sid.next = NULL;
    bsp = BioseqFind (&sid);
    if (bsp) {
      sep = SeqMgrGetSeqEntryForData (bsp);
    }
  }
  if (sep) {
    if (is_new != NULL) {
      *is_new = FALSE;
    }
    return sep;
  }
/********************/
  /*
  EntrezInit (program, TRUE, NULL);
  if (uid==0 && ! StringHasNoText (accession)) {
    if (is_na) {
      seqtype = TYP_NT;
    } else {
      seqtype = TYP_AA;
    }
    sprintf (str, "\"%s\" [ACCN]", accession);
    lsp = EntrezTLEvalString (str, seqtype, -1, NULL, NULL);
    if (lsp != NULL) {
      uid = lsp->uids [0];
    }
    LinkSetFree (lsp);
  }
  if (uid > 0) {
    sid.choice = SEQID_GI;
    sid.data.intvalue = uid;
    sid.next = NULL;
    bsp = BioseqLockById (&sid);
    if (bsp != NULL) {
      sep = SeqMgrGetSeqEntryForData (bsp);
      BioseqUnlock (bsp);
    } else {
      sep = EntrezSeqEntryGet (uid, 0);
      if (is_new != NULL) {
        *is_new = TRUE;
      }
    }
  }
  EntrezFini ();
  return sep;
  */
  return NULL;
}

static void SequinStdEditorFormMessage (ForM f, Int2 mssg)

{
  BaseFormPtr  bfp;

  bfp = (BaseFormPtr) GetObjectExtra (f);
  if (bfp != NULL) {
    switch (mssg) {
      case VIB_MSG_CLOSE :
        Remove (f);
        break;
      case VIB_MSG_QUIT :
        QuitProc ();
        break;
      case VIB_MSG_CUT :
        StdCutTextProc (NULL);
        break;
      case VIB_MSG_COPY :
        StdCopyTextProc (NULL);
        break;
      case VIB_MSG_PASTE :
        StdPasteTextProc (NULL);
        break;
      case VIB_MSG_DELETE :
        StdDeleteTextProc (NULL);
        break;
      default :
        break;
    }
  }
}

static void ProcessHelpMessage (CharPtr heading, CharPtr section)

{
  SendHelpScrollMessage (helpForm, heading, section);
}

/*
static void MakeTermListForm (void)

{
  if (! EntrezIsInited ()) {
    EntrezBioseqFetchEnable ("Sequin", TRUE);
    SequinEntrezInit ("Sequin", FALSE, NULL);
  }
  if (termListForm != NULL) return;
  termListForm = CreateTermListForm (-50, -33, "Query",
                                     TermSelectionActivateProc,
                                     TermSelectionFormMessage);
}

static void MakeDocSumForm (void)

{
  if (! EntrezIsInited ()) {
    EntrezBioseqFetchEnable ("Sequin", TRUE);
    SequinEntrezInit ("Sequin", FALSE, NULL);
  }
  if (docSumForm != NULL) return;
  docSumForm = CreateDocSumForm (-10, -90, "Document",
                                 DocumentSummaryActivateProc,
                                 DocumentSummaryFormMessage);
  if (indexerVersion) {
    UseDelayedNeighbor (docSumForm, TRUE);
  }
}
*/

void EntrezQueryProc (IteM i)

{
  /*
  MakeTermListForm ();
  MakeDocSumForm ();
  Show (termListForm);
  Select (termListForm);
  Update ();
  */
}

#if defined(OS_MAC) || defined(OS_UNIX_DARWIN)
extern Boolean Nlm_LaunchAppEx (CharPtr fileName, VoidPtr serialNumPtr, CharPtr sig);
#endif

void Entrez2QueryProc (IteM i)

{
  WatchCursor ();
  Update ();

#if defined(OS_MAC) || defined(OS_UNIX_DARWIN)
  Nlm_LaunchAppEx (NULL, NULL, "ENTZ");
#else
#if defined (OS_UNIX) || defined(WIN_MOTIF)
  system ("entrez2 &");
#else
#ifdef OS_MSWIN
  Nlm_MSWin_OpenApplication ("entrez2.exe", NULL);
#endif
#endif
#endif

  ArrowCursor ();
}

/*
*  The following callbacks process requests between forms.  This can be extended
*  to have the application track multiple term list and docsum windows.
*/

/*
static void DoRetrieveDocuments (ForM f, Int2 num, Int2 parents, Int4Ptr uids, Int2 db)

{
  if (useEntrez) {
    MakeTermListForm ();
    MakeDocSumForm ();
  }
  RetrieveDocuments (docSumForm, num, parents, uids, db);
}

static void DoRetrieveProject (ForM f, Pointer proj)

{
  if (useEntrez) {
    MakeTermListForm ();
    MakeDocSumForm ();
  }
  PointerToForm (docSumForm, proj);
}

static void DoRetrieveSimple (ForM f, ValNodePtr simpleSeqs)

{
  if (useEntrez) {
    MakeTermListForm ();
    MakeDocSumForm ();
  }
  RetrieveSimpleSeqs (docSumForm, simpleSeqs);
}

static void DoLoadNamedUidList (ForM f, CharPtr term, Int4 num, Int4Ptr uids, Int2 db)

{
  if (useEntrez) {
    MakeTermListForm ();
    MakeDocSumForm ();
  }
  LoadNamedUidList (termListForm, term, num, uids, db);
}

static void DoLaunchRecordViewer (ForM f, Int4 uid, Int2 numAlign, Int4Ptr alignuids, Int2 db)

{
  if (useEntrez) {
    MakeTermListForm ();
    MakeDocSumForm ();
  }
  LaunchRecordViewer (uid, numAlign, alignuids, db);
}

static GrouP DoMakeMedViewerLinkControls (GrouP prnt, BaseFormPtr bfp, Int2 doctype, Int4 uid)

{
  if (useEntrez) {
    if (! EntrezIsInited ()) {
      EntrezBioseqFetchEnable ("Sequin", TRUE);
      SequinEntrezInit ("Sequin", FALSE, NULL);
    }
    return MakeViewerLinkControls (prnt, bfp, doctype, uid, FALSE);
  }
  return NULL;
}

static GrouP DoMakeSeqViewerLinkControls (GrouP prnt, BaseFormPtr bfp, Int2 doctype, Int4 uid)

{
  if (useEntrez) {
    if (! EntrezIsInited ()) {
      EntrezBioseqFetchEnable ("Sequin", TRUE);
      SequinEntrezInit ("Sequin", FALSE, NULL);
    }
    return MakeViewerLinkControls (prnt, bfp, doctype, uid, TRUE);
  }
  return NULL;
}
*/

typedef struct aligngroup {
  DIALOG_MESSAGE_BLOCK
  BaseFormPtr        bfp;
  ButtoN             retrieve;
  PopuP              target;
  Int2               targetDb;
  Uint2              align_type;
  ButtoN             onlyFromThis;
} AlignGroup, PNTR AlignGroupPtr;

static void RetrieveAligns (ButtoN b)

{
  AlignGroupPtr     agp;
  BaseFormPtr       bfp;
  BioseqPtr         bsp;
  EntrezGlobalsPtr  egp;
  ValNodePtr        head;
  Int2              i;
  Int4              num;
  SeqEntryPtr       sep;
  Int4Ptr           uids;
  ValNodePtr        vnp;

  agp = (AlignGroupPtr) GetObjectExtra (b);
  if (agp == NULL) return;
  bfp = agp->bfp;
  if (bfp == NULL) return;
  egp = (EntrezGlobalsPtr) GetAppProperty ("EntrezGlobals");
  if (egp == NULL || egp->retrieveDocsProc == NULL) return;
  /*
  if (! EntrezIsInited ()) {
    SequinEntrezInit ("Sequin", FALSE, NULL);
  }
  */
  uids = NULL;
  sep = NULL;
  if (GetStatus (agp->onlyFromThis)) {
    bsp =  GetBioseqGivenIDs (bfp->input_entityID, bfp->input_itemID, bfp->input_itemtype);
    if (bsp != NULL) {
      sep = SeqMgrGetSeqEntryForData (bsp);
    }
  } else {
    sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  }
  head = GetUidsForSeqEntryAligns (sep);
  num = 0;
  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == agp->align_type) {
      num++;
    }
  }
  if (num > 0) {
    uids = MemNew ((size_t) (num + 1) * sizeof (DocUid));
    if (uids != NULL) {
      i = 0;
      for (vnp = head; i < num && vnp != NULL; vnp = vnp->next) {
        if (vnp->choice == agp->align_type) {
          uids [i] = vnp->data.intvalue;
          i++;
        }
      }
    }
    egp->retrieveDocsProc (bfp->form, num, 0, uids, agp->targetDb);
    MemFree (uids);
  }
  ValNodeFree (head);
}

static Boolean DoUpdateFetchCounts (GrouP g, SeqEntryPtr sep)

{
  AlignGroupPtr  agp;
  BaseFormPtr    bfp;
  BioseqPtr      bsp;
  ValNodePtr     head;
  Int4           num;
  Int4           rsult;
  Char           tmp [32];
  Int2           val;
  ValNodePtr     vnp;

  agp = (AlignGroupPtr) GetObjectExtra (g);
  if (agp == NULL) return FALSE;
  bfp = agp->bfp;
  if (bfp == NULL) return FALSE;
  val = GetValue (agp->target);
  num = 0;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  head = GetIdStringsForSeqEntryAligns (sep);
  rsult = ValNodeLen (head);
  if (GetStatus (agp->onlyFromThis)) {
    head = ValNodeFreeData (head);
    bsp =  GetBioseqGivenIDs (bfp->input_entityID, bfp->input_itemID, bfp->input_itemtype);
    if (bsp != NULL) {
      sep = SeqMgrGetSeqEntryForData (bsp);
      head = GetIdStringsForSeqEntryAligns (sep);
    }
  }
  agp->align_type = 0;
  agp->targetDb = TYP_NT;
  val = GetValue (agp->target);
  switch (val) {
    case 1 :
      agp->align_type = ALIGN_BLASTN;
      agp->targetDb = TYP_NT;
      break;
    case 2 :
      agp->align_type = ALIGN_BLASTP;
      agp->targetDb = TYP_AA;
      break;
    case 3 :
      agp->align_type = ALIGN_BLASTX;
      agp->targetDb = TYP_AA;
      break;
    case 4 :
      agp->align_type = ALIGN_TBLASTN;
      agp->targetDb = TYP_NT;
      break;
    default :
      break;
  }
  num = 0;
  for (vnp = head; vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == agp->align_type) {
      num++;
    }
  }
  ValNodeFreeData (head);
  sprintf (tmp, "Retrieve %ld", (long) num);
  SafeSetTitle (agp->retrieve, tmp);
  if (num > 0) {
    SafeEnable (agp->retrieve);
  } else {
    SafeDisable (agp->retrieve);
  }
  return (Boolean) (rsult > 0);
}

static void ChangeAlignTarget (PopuP p)

{
  AlignGroupPtr  agp;
  BaseFormPtr    bfp;
  SeqEntryPtr    sep;

  agp = (AlignGroupPtr) GetObjectExtra (p);
  if (agp == NULL) return;
  bfp = agp->bfp;
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  DoUpdateFetchCounts ((GrouP) agp->dialog, sep);
}

static GrouP DoMakeViewerAlignBtn (GrouP prnt, BaseFormPtr bfp)

{
  AlignGroupPtr     agp;
  EntrezGlobalsPtr  egp;
  GrouP             g;
  Boolean           macLike;
  PrompT            ppt;

  if (useEntrez) {
    /*
    if (! EntrezIsInited ()) {
      SequinEntrezInit ("Sequin", FALSE, NULL);
    }
    */
    egp = (EntrezGlobalsPtr) GetAppProperty ("EntrezGlobals");
    if (egp == NULL) return NULL;
    macLike = egp->popdownBehavior;

    agp = (AlignGroupPtr) MemNew (sizeof (AlignGroup));
    if (agp == NULL) return NULL;

    g = HiddenGroup (prnt, 5, 0, NULL);
    SetGroupSpacing (g, 10, 10);
    SetObjectExtra (g, agp, StdCleanupExtraProc);

    agp->dialog = (DialoG) g;
    agp->bfp = bfp;
    agp->retrieve = PushButton (g, "Retrieve 00000", RetrieveAligns);
    SetObjectExtra (agp->retrieve, agp, NULL);
    SetTitle (agp->retrieve, "Retrieve 0");
    SafeDisable (agp->retrieve);

    ppt = StaticPrompt (g, "Alignment:", 0, popupMenuHeight, programFont, 'l');
    agp->target = PopupList (g, macLike, (PupActnProc) ChangeAlignTarget);
    SetObjectExtra (agp->target, agp, NULL);
    PopupItem (agp->target, "BLASTN");
    PopupItem (agp->target, "BLASTP");
    PopupItem (agp->target, "BLASTX");
    PopupItem (agp->target, "TBLASTN");
    SetValue (agp->target, 1);

    agp->onlyFromThis = CheckBox (g, "Just from this sequence", (BtnActnProc) ChangeAlignTarget);
    SetObjectExtra (agp->onlyFromThis, agp, NULL);
    SetStatus (agp->onlyFromThis, TRUE);
    SafeHide (agp->onlyFromThis);

    AlignObjects (ALIGN_MIDDLE, (HANDLE) agp->retrieve, (HANDLE) ppt,
                  (HANDLE) agp->target, (HANDLE) agp->onlyFromThis, NULL);

    return g;
  }
  return NULL;
}

extern void SetupBioseqPageList (void)

{
  Char  str [32];

  seqviewprocs.pageSpecs = BioseqPageListFree (seqviewprocs.pageSpecs);
  AddBioseqPageToList (&(seqviewprocs.pageSpecs), &mapPageData);
  /* AddBioseqPageToList (&(seqviewprocs.pageSpecs), &sumPageData); */
  AddBioseqPageToList (&(seqviewprocs.pageSpecs), &asn2gphGphPageData);
  if (useOldGraphicView) {
    AddBioseqPageToList (&(seqviewprocs.pageSpecs), &gphPageData);
  }
  if (useOldAlignmentView) {
    AddBioseqPageToList (&(seqviewprocs.pageSpecs), &alnPageData);
  }
  /*
  if (useUdv) {
    AddBioseqPageToList (&(seqviewprocs.pageSpecs), &udvPageData);
    AddBioseqPageToList (&(seqviewprocs.pageSpecs), &ddvPageData);
  } else {
    AddBioseqPageToList (&(seqviewprocs.pageSpecs), &seqPageData);
  }
  */
  AddBioseqPageToList (&(seqviewprocs.pageSpecs), &seqpnlPageData);
  AddBioseqPageToList (&(seqviewprocs.pageSpecs), &seqAlnPnlPageData);
  if (useOldSequenceView) {
    AddBioseqPageToList (&(seqviewprocs.pageSpecs), &seqPageData);
  }
  AddBioseqPageToList (&(seqviewprocs.pageSpecs), &gbgnPageData);
  if (Nlm_GetAppProperty ("SequinUseEMBLStyle") != NULL) {
    AddBioseqPageToList (&(seqviewprocs.pageSpecs), &emblPageData);
    AddBioseqPageToList (&(seqviewprocs.pageSpecs), &gnbkPageData);
    if (GetSequinAppParam ("SETTINGS", "NUCPAGE", "EMBL", str, sizeof (str))) {
      seqviewprocs.initNucLabel = MemFree (seqviewprocs.initNucLabel);
      seqviewprocs.initNucLabel = StringSaveNoNull (str);
    }
  } else if (Nlm_GetAppProperty ("SequinUseDDBJStyle") != NULL) {
    AddBioseqPageToList (&(seqviewprocs.pageSpecs), &ddbjPageData);
    AddBioseqPageToList (&(seqviewprocs.pageSpecs), &emblPageData);
    if (GetSequinAppParam ("SETTINGS", "NUCPAGE", "DDBJ", str, sizeof (str))) {
      seqviewprocs.initNucLabel = MemFree (seqviewprocs.initNucLabel);
      seqviewprocs.initNucLabel = StringSaveNoNull (str);
    }
  } else {
    AddBioseqPageToList (&(seqviewprocs.pageSpecs), &gnbkPageData);
    AddBioseqPageToList (&(seqviewprocs.pageSpecs), &emblPageData);
    if (GetSequinAppParam ("SETTINGS", "NUCPAGE", "GenBank", str, sizeof (str))) {
      seqviewprocs.initNucLabel = MemFree (seqviewprocs.initNucLabel);
      seqviewprocs.initNucLabel = StringSaveNoNull (str);
    }
  }
  AddBioseqPageToList (&(seqviewprocs.pageSpecs), &gnptPageData);
  AddBioseqPageToList (&(seqviewprocs.pageSpecs), &ftblPageData);
  AddBioseqPageToList (&(seqviewprocs.pageSpecs), &fstaPageData);
  AddBioseqPageToList (&(seqviewprocs.pageSpecs), &qualPageData);
  AddBioseqPageToList (&(seqviewprocs.pageSpecs), &asnPageData);
  AddBioseqPageToList (&(seqviewprocs.pageSpecs), &xmlPageData);
  AddBioseqPageToList (&(seqviewprocs.pageSpecs), &gbseqPageData);

  AddBioseqPageToList (&(seqviewprocs.pageSpecs), &dskPageData);
}

static Boolean DeltaLitOnly (BioseqPtr bsp)

{
  ValNodePtr  vnp;

  if (bsp == NULL || bsp->repr != Seq_repr_delta) return FALSE;
  for (vnp = (ValNodePtr)(bsp->seq_ext); vnp != NULL; vnp = vnp->next) {
    if (vnp->choice == 1) return FALSE;
  }
  return TRUE;
}

static Int2 LIBCALLBACK DeltaToRawConvertFunc (Pointer data)

{
  BioseqPtr         bsp = NULL;
  OMProcControlPtr  ompcp;
  
  if (indexerVersion)
  {
    SeqEditFunc (data);
  }
  else
  {
    ompcp = (OMProcControlPtr) data;
    if (ompcp == NULL || ompcp->proc == NULL) 
      return OM_MSG_RET_ERROR;
    switch (ompcp->input_itemtype) {
      case OBJ_BIOSEQ :
        bsp = (BioseqPtr) ompcp->input_data;
        break;
      case 0 :
        return OM_MSG_RET_ERROR;
      default :
        return OM_MSG_RET_ERROR;
    }
    if (bsp == NULL) {
       return OM_MSG_RET_ERROR;
    }
    if (! DeltaLitOnly (bsp)) {
      return OM_MSG_RET_OK;
    }
    if (Message (MSG_YN, "Convert near delta Bioseq to raw Bioseq?") == ANS_NO) {
      return OM_MSG_RET_DONE;
    }
    SegOrDeltaBioseqToRaw (bsp);
    Message (MSG_OK, "Converted to raw, now launch editor again");
  }
  return OM_MSG_RET_DONE;
}

#define REGISTER_DELTA_BIOSEQ_EDIT ObjMgrProcLoad(OMPROC_EDIT,"Delta Bioseq Convert","DeltaBioseqConverter",OBJ_BIOSEQ,Seq_repr_delta,OBJ_BIOSEQ,Seq_repr_raw,NULL,DeltaToRawConvertFunc,PROC_PRIORITY_DEFAULT)

static void SetupDesktop (void)

{
  Boolean     allowalign;
  Boolean     cdsMrnaOneToOne;
  FeatDefPtr  curr;
  Uint1       key;
  CharPtr     label = NULL;
  Char        proclabel [64];
  Char        procname [64];
  Boolean     readOnlyDbTags;
  Char        str [PATH_MAX];
  Uint2       subtype;
  Int2        val;
  Boolean     validateExons;

  SetAppProperty ("SequinAppVersion", SEQUIN_APPLICATION);

  MemSet ((Pointer) (&medviewprocs), 0, sizeof (MedlineViewProcs));
  medviewprocs.cleanupObjectPtr = FALSE;
  medviewprocs.activateForm = MedlineViewFormActivated;
  medviewprocs.closeForm = NULL;
  medviewprocs.useFolderTabs = CHANGE_VIEW_POPUP;
  /*
  medviewprocs.useFolderTabs = CHANGE_VIEW_RADIOBUTTONS;
  medviewprocs.initPage = CITATION_PAGE;
  */
#ifndef WIN_MAC
  medviewprocs.createMenus = MedlineViewFormMenus;
#endif
  medviewprocs.showAsnPage = TRUE;
  if (indexerVersion) {
    medviewprocs.useScrollText = TRUE;
  } else {
    medviewprocs.useScrollText = FALSE;
  }
  medviewprocs.handleMessages = SequinMedlineFormMessage;
  /*
  medviewprocs.makeControls = DoMakeMedViewerLinkControls;
  */
  SetAppProperty ("MedlineDisplayForm", &medviewprocs);

  MemSet ((Pointer) (&seqviewprocs), 0, sizeof (SeqViewProcs));
  seqviewprocs.hasTargetControl = TRUE;
  seqviewprocs.hasDoneButton = TRUE;
  seqviewprocs.launchEditors = TRUE;
  seqviewprocs.launchSubviewers = FALSE;
  seqviewprocs.sendSelectMessages = TRUE;
  seqviewprocs.highlightSelections = TRUE;
  seqviewprocs.cleanupObjectPtr = FALSE;
  seqviewprocs.activateForm = BioseqViewFormActivated;
  seqviewprocs.closeForm = StdSendCloseWindowMessageProc;
  seqviewprocs.useFolderTabs = CHANGE_VIEW_POPUP;
  /*
  seqviewprocs.initNucPage = NUCASN2FF_PAGE_1;
  seqviewprocs.initProtPage = PROTGENPEPT_PAGE;
  */
#ifndef WIN_MAC
  seqviewprocs.createMenus = BioseqViewFormMenus;
#endif
#ifdef WIN_MOTIF
  if (indexerVersion) {
    if (GetSequinAppParam ("SETTINGS", "WGS", NULL, str, sizeof (str))
        && StringICmp (str, "TRUE") == 0) {
      seqviewprocs.createToolBar = BioseqViewFormWGSToolBar;
    } else {
      seqviewprocs.createToolBar = BioseqViewFormToolBar;
    }
  }
#endif
#ifdef WIN_MSWIN
  if (indexerVersion) {
    if (GetSequinAppParam ("SETTINGS", "WGS", NULL, str, sizeof (str))
        && StringICmp (str, "TRUE") == 0) {
      seqviewprocs.createToolBar = BioseqViewFormWGSToolBar;
    } else {
      seqviewprocs.createToolBar = BioseqViewFormToolBar;
    }
  }
#endif
/*#ifdef INTERNAL_NCBI_SEQUIN*/
  if (indexerVersion) {
    seqviewprocs.allowScrollText = TRUE;
    seqviewprocs.startInScrollText = FALSE;
  } else {
    seqviewprocs.allowScrollText = FALSE;
    seqviewprocs.startInScrollText = FALSE;
  }
/*#endif*/
  seqviewprocs.handleMessages = SequinSeqViewFormMessage;
  /*
  seqviewprocs.makeControls = DoMakeSeqViewerLinkControls;
  seqviewprocs.updateControls = UpdateViewerLinkTarget;
  */
  seqviewprocs.makeAlignBtn = DoMakeViewerAlignBtn;
  seqviewprocs.updateCounts = DoUpdateFetchCounts;
  seqviewprocs.filepath = NULL;

  SetupBioseqPageList ();

  SetAppProperty ("SeqDisplayForm", &seqviewprocs);

  MemSet ((Pointer) (&seqedprocs), 0, sizeof (SeqEditViewProcs));
#ifdef WIN_MAC
  seqedprocs.activateForm = SeqEditFormActivated;
#endif
  seqedprocs.handleMessages = SequinSeqEditFormMessage;
  seqedprocs.minPixelWidth = 0;
  seqedprocs.minPixelHeight = 0;
  seqedprocs.showfeat = FALSE;
  seqedprocs.extended_align_menu = FALSE;
  seqedprocs.extended_dist_menu = FALSE;
  seqedprocs.extended_tree_menu = FALSE;

  if (useEntrez) {
    seqedprocs.download = SeqEdDownload;
  }

  SetAppProperty ("SeqEditDisplayForm", &seqedprocs);

  MemSet ((Pointer) (&stdedprocs), 0, sizeof (StdEditorProcs));
#ifdef WIN_MAC
  stdedprocs.activateForm = StdEditorFormActivated;
#endif
  stdedprocs.handleMessages = SequinStdEditorFormMessage;
  stdedprocs.duplicateExisting = FALSE;
  SetAppProperty ("StdEditorForm", &stdedprocs);

  MemSet ((Pointer) (&valdtrprocs), 0, sizeof (StdEditorProcs));
#ifdef WIN_MAC
  valdtrprocs.activateForm = StdValidatorFormActivated;
#endif
  valdtrprocs.handleMessages = SequinStdEditorFormMessage;
  valdtrprocs.duplicateExisting = FALSE;
  SetAppProperty ("StdValidatorForm", &valdtrprocs);

  MemSet ((Pointer) (&txtviewprocs), 0, sizeof (TextViewProcs));
#ifdef WIN_MAC
  txtviewprocs.activateForm = TextViewProcFormActivated;
#endif
  if (indexerVersion) {
    txtviewprocs.useScrollText = TRUE;
  } else {
    txtviewprocs.useScrollText = FALSE;
  }
  SetAppProperty ("TextDisplayForm", &txtviewprocs);

  SetAppProperty ("HelpMessageProc", (Pointer) ProcessHelpMessage);

  MemSet ((Pointer) (&pubedprocs), 0, sizeof (PubdescEditProcs));
  if (newMedarch) {
    pubedprocs.lookupArticle = LookupAnArticleFuncNew;
    pubedprocs.lookupJournal = LookupJournalFuncNew;
  } else if (useMedarch) {
    pubedprocs.lookupArticle = LookupAnArticleFunc;
    pubedprocs.lookupJournal = LookupJournalFunc;
  }
/*#ifdef REPLACE_THIS*/
  if (indexerVersion) {
    pubedprocs.replaceThis = TRUE;
  }
/*#endif*/
  SetAppProperty ("PubdescEditForm", &pubedprocs);

  MemSet ((Pointer) (&biosrcedprocs), 0, sizeof (BioSourceEditProcs));
  if (useTaxon) {
    biosrcedprocs.lookupTaxonomy = LookupTaxonomyFunc;
  }
  SetAppProperty ("BioSourcEditForm", &biosrcedprocs);

  readOnlyDbTags = TRUE;

/*#ifdef INTERNAL_NCBI_SEQUIN*/
  if (indexerVersion) {
    SetAppProperty ("InternalNcbiSequin", (void *) 1024);
    readOnlyDbTags = FALSE;
  }
/*#endif*/

  if (genomeCenter != NULL) {
    SetAppProperty ("GenomeCenterSequin", (void *) 1024);
    readOnlyDbTags = FALSE;
  }

  if (readOnlyDbTags) {
    SetAppProperty ("ReadOnlyDbTags", (void *) 1024);
  }

  SetAppProperty ("NewSequinGraphicalViewer", (void *) 1024);
  SetAppProperty ("NewSequinLayoutOverride", (void *) 1024);

  if (GetSequinAppParam ("SETTINGS", "BROWSER", NULL, str, sizeof (str))) {
    SetAppProperty ("MedviewBrowserPath", (void *) StringSave (str));
  }

  if (GetSequinAppParam ("PREFERENCES", "MINPIXELWIDTH", NULL, str, sizeof (str))) {
    if (StrToInt (str, &val) && val > 0) {
      val = MIN (val, screenRect.right);
      medviewprocs.minPixelWidth = val;
      seqviewprocs.minPixelWidth = val;
    }
  }

  if (GetSequinAppParam ("PREFERENCES", "MINPIXELHEIGHT", NULL, str, sizeof (str))) {
    if (StrToInt (str, &val) && val > 0) {
      val = MIN (val, screenRect.bottom);
      medviewprocs.minPixelHeight = val;
      seqviewprocs.minPixelHeight = val;
    }
  }

  if (GetSequinAppParam ("PREFERENCES", "SEQEDPIXELWIDTH", NULL, str, sizeof (str))) {
    if (StrToInt (str, &val) && val > 0) {
      val = MIN (val, screenRect.right);
      seqedprocs.minPixelWidth = val;
    }
  }

  if (GetSequinAppParam ("PREFERENCES", "SEQEDPIXELHEIGHT", NULL, str, sizeof (str))) {
    if (StrToInt (str, &val) && val > 0) {
      val = MIN (val, screenRect.bottom);
      seqedprocs.minPixelHeight = val;
    }
  }

  if (indexerVersion) {
    if (seqviewprocs.minPixelWidth < 500) {
      seqviewprocs.minPixelWidth = 500;
    }
    if (seqviewprocs.minPixelHeight < 600) {
      seqviewprocs.minPixelHeight = 600;
    }
    if (txtviewprocs.minPixelWidth < 750) {
      txtviewprocs.minPixelWidth = 750;
    }
    if (txtviewprocs.minPixelHeight < 600) {
      txtviewprocs.minPixelHeight = 600;
    }
  }

  if (screenRect.bottom < 780) {
    if (seqviewprocs.minPixelHeight > 500) {
      seqviewprocs.minPixelHeight = 500;
    }
    if (txtviewprocs.minPixelHeight > 500) {
      txtviewprocs.minPixelHeight = 500;
    }
  }

  if (GetSequinAppParam ("PREFERENCES", "TEXTPIXELWIDTH", NULL, str, sizeof (str))) {
    if (StrToInt (str, &val) && val > 0) {
      val = MIN (val, screenRect.right);
      txtviewprocs.minPixelWidth = val;
    }
  }

  if (GetSequinAppParam ("PREFERENCES", "TEXTPIXELHEIGHT", NULL, str, sizeof (str))) {
    if (StrToInt (str, &val) && val > 0) {
      val = MIN (val, screenRect.bottom);
      txtviewprocs.minPixelHeight = val;
    }
  }

  if (GetSequinAppParam ("PREFERENCES", "SEQEDSHOWFEAT", NULL, str, sizeof (str))) {
    if (StringICmp (str, "TRUE") == 0) {
      seqedprocs.showfeat = TRUE;
    }
  }

  if (GetSequinAppParam ("PREFERENCES", "SEQEDEXTENDALIGN", NULL, str, sizeof (str))) {
    if (StringICmp (str, "TRUE") == 0) {
      seqedprocs.extended_align_menu = TRUE;
    }
  }

  if (GetSequinAppParam ("PREFERENCES", "SEQEDEXTENDDIST", NULL, str, sizeof (str))) {
    if (StringICmp (str, "TRUE") == 0) {
      seqedprocs.extended_dist_menu = TRUE;
    }
  }

  if (GetSequinAppParam ("PREFERENCES", "SEQEDEXTENDTREE", NULL, str, sizeof (str))) {
    if (StringICmp (str, "TRUE") == 0) {
      seqedprocs.extended_tree_menu = TRUE;
    }
  }

  if (GetSequinAppParam ("PREFERENCES", "USEOLDASN", NULL, str, sizeof (str))) {
    if (StringICmp (str, "TRUE") == 0) {
      leaveAsOldAsn = TRUE;
    }
  }

  if (GetSequinAppParam ("SETTINGS", "XMLPREFIX", NULL, str, sizeof (str))) {
    AsnSetXMLmodulePrefix (StringSave (str));
  }

  validateExons = TRUE;
  if (GetSequinAppParam ("PREFERENCES", "VALIDATEEXONS", NULL, str, sizeof (str))) {
    if (StringICmp (str, "FALSE") == 0) {
      validateExons = FALSE;
    }
  }
  if (validateExons) {
    SetAppProperty ("ValidateExons", (void *) 1024);
  }
  /*
  SetAppProperty ("SpliceValidateAsError", (void *) 1024);
  */

  cdsMrnaOneToOne = FALSE;
  if (GetSequinAppParam ("PREFERENCES", "VALIDATECDSMRNA", NULL, str, sizeof (str))) {
    if (StringICmp (str, "TRUE") == 0) {
      cdsMrnaOneToOne = FALSE;
    }
  }
  if (cdsMrnaOneToOne) {
    SetAppProperty ("ValidateCDSmRNAoneToOne", (void *) 1024);
  }

  testLatLonSubregion = FALSE;
  if (GetSequinAppParam ("PREFERENCES", "VALIDATELATLONSUBREGION", NULL, str, sizeof (str))) {
    if (StringICmp (str, "TRUE") == 0) {
      testLatLonSubregion = TRUE;
    }
  }

  strictLatLonCountry = FALSE;
  if (GetSequinAppParam ("PREFERENCES", "VALIDATELATLONSTRICT", NULL, str, sizeof (str))) {
    if (StringICmp (str, "TRUE") == 0) {
      strictLatLonCountry = TRUE;
    }
  }

  /*
  if (GetSequinAppParam ("SETTINGS", "NUCFIELD", NULL, str, sizeof (str))) {
    if (StrToInt (str, &val) && val >= 0) {
      seqviewprocs.initNucPage = val;
    }
  }

  if (GetSequinAppParam ("SETTINGS", "PRTFIELD", NULL, str, sizeof (str))) {
    if (StrToInt (str, &val) && val >= 0) {
      seqviewprocs.initProtPage = val;
    }
  }
  */

  if (GetSequinAppParam ("SETTINGS", "TARGETCONTROL", NULL, str, sizeof (str))) {
    if (StrToInt (str, &val) && val >= 0) {
      seqviewprocs.useFolderTabs = val;
    }
  }

  if (GetSequinAppParam ("SETTINGS", "ALLOWSCROLLTEXT", NULL, str, sizeof (str))) {
    if (StringICmp (str, "TRUE") == 0) {
      seqviewprocs.allowScrollText = TRUE;
      seqviewprocs.startInScrollText = FALSE;
    }
  }

  if (GetSequinAppParam ("SETTINGS", "MEDPAGE", "Abstract", str, sizeof (str))) {
    medviewprocs.initMedLabel = StringSaveNoNull (str);
  }
  if (GetSequinAppParam ("SETTINGS", "NUCPAGE", "GenBank", str, sizeof (str))) {
    seqviewprocs.initNucLabel = StringSaveNoNull (str);
  }
  if (GetSequinAppParam ("SETTINGS", "PRTPAGE", "GenPept", str, sizeof (str))) {
    seqviewprocs.initProtLabel = StringSaveNoNull (str);
  }
  if (GetSequinAppParam ("SETTINGS", "GENMPAGE", "Map", str, sizeof (str))) {
    seqviewprocs.initGenomeLabel = StringSaveNoNull (str);
  }

  if (useEntrez) {
    seqviewprocs.lockFarComponents = TRUE;
  }
  if (GetSequinAppParam ("SETTINGS", "LOCKFAR", NULL, str, sizeof (str))) {
    if (StringICmp (str, "FALSE") == 0) {
      seqviewprocs.lockFarComponents = FALSE;
    }
  }

  if (GetSequinAppParam ("PREFERENCES", "GPHVIEWSCOREALIGNS", NULL, str, sizeof (str))) {
    if (StringICmp (str, "TRUE") == 0) {
      gphviewscorealigns = TRUE;
    }
  }
  if (gphviewscorealigns) {
    SetAppProperty("GPHVIEWSCOREALIGNS", (void *) 1);
  }

  MemSet ((Pointer) (&entrezglobals), 0, sizeof (EntrezGlobals));
  /*
  entrezglobals.retrieveDocsProc = DoRetrieveDocuments;
  entrezglobals.retrieveProjectProc = DoRetrieveProject;
  entrezglobals.retrieveSimpleProc = DoRetrieveSimple;
  entrezglobals.loadNamedUidProc = DoLoadNamedUidList;
  entrezglobals.launchViewerProc = DoLaunchRecordViewer;
#ifndef WIN_MAC
  entrezglobals.createTrmLstMenus = TermListFormMenus;
  entrezglobals.createDocSumMenus = DocSumFormMenus;
#endif
  SetAppProperty ("EntrezGlobals", &entrezglobals);
  */

  entrezglobals.showAsn = TRUE;

  entrezglobals.persistDefault = TRUE;
  if (GetSequinAppParam ("PREFERENCES", "PARENTSPERSIST", NULL, str, sizeof (str) - 1)) {
    if (StringICmp (str, "FALSE") == 0) {
      entrezglobals.persistDefault = FALSE;
    }
  }
  entrezglobals.alignDefault = TRUE;
  if (GetSequinAppParam ("PREFERENCES", "ALIGNCHECKED", NULL, str, sizeof (str) - 1)) {
    if (StringICmp (str, "FALSE") == 0) {
      entrezglobals.alignDefault = FALSE;
    }
  }
  seqviewprocs.alignDefault = entrezglobals.alignDefault;
  if (GetSequinAppParam ("PREFERENCES", "POPDOWN", NULL, str, sizeof (str) - 1)) {
    if (StringICmp (str, "TRUE") == 0) {
      entrezglobals.popdownBehavior = TRUE;
    }
  }
  if (GetSequinAppParam ("PREFERENCES", "LOOKUPDIRECT", NULL, str, sizeof (str) - 1)) {
    if (StringICmp (str, "TRUE") == 0) {
      entrezglobals.lookupDirect = TRUE;
    }
  }
  entrezglobals.sortFields = TRUE;
  if (GetSequinAppParam ("PREFERENCES", "SORTFIELDS", NULL, str, sizeof (str) - 1)) {
    if (StringICmp (str, "FALSE") == 0) {
      entrezglobals.sortFields = FALSE;
    }
  }
  if (GetSequinAppParam ("SETTINGS", "DATABASE", "MEDLINE", str, sizeof (str))) {
    entrezglobals.initDatabase = StringSaveNoNull (str);
  }
  if (GetSequinAppParam ("SETTINGS", "FIELD", "All Fields", str, sizeof (str))) {
    entrezglobals.initField = StringSaveNoNull (str);
  }
  if (GetSequinAppParam ("SETTINGS", "MODE", "Automatic", str, sizeof (str))) {
    entrezglobals.initMode = StringSaveNoNull (str);
  }

/*#ifdef INTERNAL_NCBI_SEQUIN*/
  if (indexerVersion) {
    allowalign = TRUE;
  } else {
/*#else*/
    allowalign = FALSE;
  }
/*#endif*/
  if (GetSequinAppParam ("SETTINGS", "ALLOWALIGNMENT", NULL, str, sizeof (str))) {
    if (StringICmp (str, "TRUE") == 0) {
      allowalign = TRUE;
    }
  }
  if (allowalign) {
    SetAppProperty ("AllowAlignment", (void *) 1024);
  }

  /* register types and functions first */

  REGISTER_MEDLINE_VIEW;
  if (smartnetMode) {
    REGISTER_SMART_SEQENTRY_VIEW;
  } else {
    REGISTER_NEW_SEQENTRY_VIEW;
  }

  curr = FeatDefFindNext (NULL, &key, &label, FEATDEF_ANY, TRUE);
  while (curr != NULL) {
    if (key != FEATDEF_BAD && curr->seqfeat_key == SEQFEAT_IMP) {
      subtype = curr->featdef_key;
      if (subtype != FEATDEF_source) {
        sprintf (procname, "Edit %s", curr->menulabel);
        StringNCpy_0 (proclabel, curr->typelabel, sizeof (proclabel));
        if (proclabel [0] == '-') {
          proclabel [0] = '~';
        }
        REGISTER_IMPORT_EDIT(procname,proclabel,subtype);
      }
    }
    curr = FeatDefFindNext (curr, &key, &label, FEATDEF_ANY, TRUE);
  }

  curr = FeatDefFindNext (NULL, &key, &label, FEATDEF_ANY, TRUE);
  while (curr != NULL) {
    if (key != FEATDEF_BAD && curr->seqfeat_key == SEQFEAT_RNA) {
      subtype = curr->featdef_key;
      sprintf (procname, "Edit %s", curr->menulabel);
      StringNCpy_0 (proclabel, curr->typelabel, sizeof (proclabel));
      if (proclabel [0] == '-') {
        proclabel [0] = '~';
      }
      REGISTER_RNA_EDIT(procname,proclabel,subtype);
    }
    curr = FeatDefFindNext (curr, &key, &label, FEATDEF_ANY, TRUE);
  }

  curr = FeatDefFindNext (NULL, &key, &label, FEATDEF_ANY, TRUE);
  while (curr != NULL) {
    if (key != FEATDEF_BAD && curr->seqfeat_key == SEQFEAT_PROT) {
      subtype = curr->featdef_key;
      sprintf (procname, "Edit %s", curr->menulabel);
      StringNCpy_0 (proclabel, curr->typelabel, sizeof (proclabel));
      if (proclabel [0] == '-') {
        proclabel [0] = '~';
      }
      REGISTER_PROT_EDIT(procname,proclabel,subtype);
    }
    curr = FeatDefFindNext (curr, &key, &label, FEATDEF_ANY, TRUE);
  }

  REGISTER_GENE_EDIT;
  REGISTER_CDRGN_EDIT;
  REGISTER_RGNFEAT_EDIT;
  REGISTER_CMNTFEAT_EDIT;
  REGISTER_BIOSOURCE_FEAT_EDIT;
  REGISTER_BIOSOURCE_DESC_EDIT;
  REGISTER_MOLINFO_EDIT;
  REGISTER_PUBDESC_FEAT_EDIT;
  REGISTER_PUBDESC_DESC_EDIT;
  REGISTER_TITLE_EDIT;
  REGISTER_COMMENT_EDIT;
  REGISTER_NAME_EDIT;
  REGISTER_REGION_EDIT;
  REGISTER_BOND_EDIT;
  REGISTER_SITE_EDIT;
  REGISTER_PSEC_EDIT;
  REGISTER_GENBANK_EDIT;
  REGISTER_CREATE_DATE_EDIT;
  REGISTER_UPDATE_DATE_EDIT;

  if (indexerVersion) {
    REGISTER_INGENUE;
  }
  REGISTER_NEW_BIOSEQ_EDIT;
  REGISTER_DELTA_BIOSEQ_EDIT;
  REGISTER_NEW_SEQALIGN_EDIT;
  REGISTER_NEW_SEQANNOT_EDIT;
  REGISTER_NEW_SEQALIGN_VIEW; 
  REGISTER_BIOSEQ_SEG_EDIT;
  REGISTER_BIOSEQ_SET_EDIT;
  REGISTER_SUBMITBLOCK_EDIT;
  REGISTER_SEQSUBCIT_EDIT;

  SetupSequinFilters ();
}

static FonT ChooseAFont (CharPtr param, CharPtr dfault)

{
  FonT  f;
  Char  str [128];

  f = NULL;
  if (GetSequinAppParam ("FONTS", param, NULL, str, sizeof (str))) {
    f = ParseFont (str);
  } else {
    if (! indexerVersion) {
      SetSequinAppParam ("FONTS", param, dfault);
    }
    f = ParseFont (dfault);
  }
  return f;
}

static void SetupCommonFonts (void)

{
#ifdef WIN_MAC
  medviewprocs.jourfnt = ChooseAFont ("JOURNAL", "Geneva,10,i");
  medviewprocs.volfnt = ChooseAFont ("VOLUME", "Geneva,10,b");
  medviewprocs.pagesfnt = ChooseAFont ("PAGES", "Geneva,10");
  medviewprocs.titlefnt = ChooseAFont ("TITLE", "Times,14,b");
  medviewprocs.authorsfnt = ChooseAFont ("AUTHORS", "Times,14");
  medviewprocs.affilfnt = ChooseAFont ("AFFILIATION", "Times,12");
  medviewprocs.abstractfnt = ChooseAFont ("ABSTRACT", "Geneva,10");
  medviewprocs.meshfnt = ChooseAFont ("MESH", "Monaco,9");
  medviewprocs.displayFont = ChooseAFont ("DISPLAY", "Monaco,9");
#endif
#ifdef WIN_MSWIN
  medviewprocs.jourfnt = ChooseAFont ("JOURNAL", "Arial,11,i");
  medviewprocs.volfnt = ChooseAFont ("VOLUME", "Arial,11,b");
  medviewprocs.pagesfnt = ChooseAFont ("PAGES", "Arial,11");
  medviewprocs.titlefnt = ChooseAFont ("TITLE", "Times New Roman,14,b");
  medviewprocs.authorsfnt = ChooseAFont ("AUTHORS", "Times New Roman,14");
  medviewprocs.affilfnt = ChooseAFont ("AFFILIATION", "Times New Roman,11");
  medviewprocs.abstractfnt = ChooseAFont ("ABSTRACT", "Times New Roman,11");
  medviewprocs.meshfnt = ChooseAFont ("MESH", "Times New Roman,9");
  medviewprocs.displayFont = ChooseAFont ("DISPLAY", "Courier New,10");
#endif
#ifdef WIN_MOTIF
  medviewprocs.jourfnt = ChooseAFont ("JOURNAL", "Helvetica,12,i");
  medviewprocs.volfnt = ChooseAFont ("VOLUME", "Helvetica,12,b");
  medviewprocs.pagesfnt = ChooseAFont ("PAGES", "Helvetica,12");
  medviewprocs.titlefnt = ChooseAFont ("TITLE", "Times,18,b");
  medviewprocs.authorsfnt = ChooseAFont ("AUTHORS", "Times,18");
  medviewprocs.affilfnt = ChooseAFont ("AFFILIATION", "Times,14");
  medviewprocs.abstractfnt = ChooseAFont ("ABSTRACT", "Times,14");
  medviewprocs.meshfnt = ChooseAFont ("MESH", "Times,12");
  medviewprocs.displayFont = ChooseAFont ("DISPLAY", "Courier,10");
#endif

#ifdef WIN_MAC
  seqviewprocs.displayFont = ChooseAFont ("DISPLAY", "Monaco,9");
#endif
#ifdef WIN_MSWIN
  if (indexerVersion) {
    seqviewprocs.displayFont = ChooseAFont ("DISPLAY", "Courier New,11");
  } else {
    seqviewprocs.displayFont = ChooseAFont ("DISPLAY", "Courier New,10");
  }
#endif
#ifdef WIN_MOTIF
  seqviewprocs.displayFont = ChooseAFont ("DISPLAY", "Courier,10");
#endif

#ifdef WIN_MAC
  entrezglobals.docsumFont = ChooseAFont ("FETCHED", "Monaco,9");
#endif
#ifdef WIN_MSWIN
  entrezglobals.docsumFont = ChooseAFont ("FETCHED", "Courier New,10");
#endif
#ifdef WIN_MOTIF
  entrezglobals.docsumFont = ChooseAFont ("FETCHED", "Courier,10");
#endif

#ifdef WIN_MAC
  txtviewprocs.displayFont = ChooseAFont ("DISPLAY", "Monaco,9");
#endif
#ifdef WIN_MSWIN
  txtviewprocs.displayFont = ChooseAFont ("DISPLAY", "Courier New,10");
#endif
#ifdef WIN_MOTIF
  txtviewprocs.displayFont = ChooseAFont ("DISPLAY", "Courier,10");
#endif
}

static void RevStringUpper (CharPtr str)
{
	CharPtr nd;
	Char tmp;

		if (str == NULL)
			return;
    nd = str;
	while (*nd != '\0')
		nd++;
	nd--;

	while (nd > str)
	{
		tmp = TO_UPPER(*nd);
		*nd = TO_UPPER(*str);
		*str = tmp;
		nd--; str++;
	}

	if (nd == str)
		*nd = TO_UPPER(*nd);
	return;
}

static void ObjMgrReport (IteM i)

{
  BioseqPtr                  bsp;
  Char                       buf [41];
  FILE                       *fp;
  Uint4                      j, num;
  ObjMgrDataPtr              omdp;
  ObjMgrDataPtr PNTR         omdpp;
  ObjMgrPtr                  omp;
  Char                       path [PATH_MAX];
  SeqIdIndexElementPtr PNTR  sipp;
  SeqMgrPtr                  smp;

  omp = ObjMgrGet ();
  if (omp == NULL) return;
  smp = SeqMgrGet ();
  if (smp == NULL) return;
  TmpNam (path);
  fp = FileOpen (path, "w");
  if (fp == NULL) return;
  fprintf (fp, "totobj %d, currobj %d\n", (int) omp->totobj, (int) omp->currobj);
  omdpp = omp->datalist;
  if (omdpp != NULL) {
    num = omp->currobj;
    for (j = 0; j < num; j++) {
      omdp = omdpp [j];
      if (omdp != NULL && omdp->datatype == OBJ_BIOSEQ) {
        bsp = (BioseqPtr) omdp->dataptr;
        if (bsp != NULL) {
          SeqIdWrite (bsp->id, buf, PRINTID_REPORT, sizeof (buf) - 1);
          fprintf (fp, "%4ld %s\n", (long) j, buf);
        }
      }
    }
  }
  fprintf (fp, "seqid index count %ld\n", (long) smp->BioseqIndexCnt);
  num = smp->BioseqIndexCnt;
  sipp = smp->BioseqIndex;
  if (sipp != NULL) {
    for (j = 0; j < num; j++) {
      StringNCpy_0 (buf, sipp [j]->str, sizeof (buf));
      RevStringUpper (buf);
      fprintf (fp, "%4ld %s\n", (long) j, buf);
    }
  }
  FileClose (fp);
  LaunchGeneralTextViewer (path, "Object Manager Report");
  FileRemove (path);
}

#ifdef WIN_MAC

static void ToggleOldAsnProc (IteM i)

{
  leaveAsOldAsn = GetStatus (i);
}

static void SetupMacMenus (void)

{
  MenU  m;
  Int2  mssgadd;
  Int2  mssgalign;
  Int2  mssgdup;
  Int2  mssgfeatprop;
  Int2  mssgseq;
  Int2  mssgsub;
  Int2  mssgupd, mssgupd_idx = 0;
  Int2  mssgext;
  /*
  MenU  sub;
  */

  m = AppleMenu (NULL);
  AddAboutAndHelpMenuItems (m);
  DeskAccGroup (m);

  m = PulldownMenu (NULL, "File");
  openItem = CommandItem (m, "Open.../O", MacReadNewAsnProc);
  SetAppProperty (SEQFORM_OPEN_ITEM, openItem);
  closeItem = FormCommandItem (m, "Close", NULL, VIB_MSG_CLOSE);
  SetAppProperty (SEQFORM_CLOSE_ITEM, closeItem);
  
  SetAppProperty (SEQFORM_INIT_ACTIVE, &initialFormsActive);
  
  SeparatorItem (m);
  importItem = FormCommandItem (m, "Import.../I", NULL, VIB_MSG_IMPORT);
  exportItem = FormCommandItem (m, "Export.../E", NULL, VIB_MSG_EXPORT);
  SeparatorItem (m);
  duplicateViewItem = CommandItem (m, "Duplicate View", DuplicateViewProc);
  SeparatorItem (m);
  saveItem = FormCommandItem (m, "Save/S", NULL, VIB_MSG_SAVE);
  saveAsItem = FormCommandItem (m, "Save As...", NULL, VIB_MSG_SAVE_AS);
  SeparatorItem (m);
  restoreItem = CommandItem (m, "Restore...", RestoreSeqEntryProc);
  SeparatorItem (m);
  prepareItem = CommandItem (m, "Prepare Submission...", PrepareSeqSubmitProc);
  /*
  submitItem = CommandItem (m, "Submit to NCBI", SubmitToNCBI);
  */
  /*
  SeparatorItem (m);
  CommandItem (m, "Propagate Top Descriptors", ForcePropagate);
  */
  SeparatorItem (m);
  printItem = FormCommandItem (m, "Print", NULL, VIB_MSG_PRINT);
  SeparatorItem (m);
  FormCommandItem (m, "Quit/Q", NULL, VIB_MSG_QUIT);

  m = PulldownMenu (NULL, "Edit");
  undoItem = FormCommandItem (m, UNDO_MENU_ITEM, NULL, VIB_MSG_UNDO);
  Disable (undoItem);
  SeparatorItem (m);
  cutItem = FormCommandItem (m, CUT_MENU_ITEM, NULL, VIB_MSG_CUT);
  copyItem = FormCommandItem (m, COPY_MENU_ITEM, NULL, VIB_MSG_COPY);
  pasteItem = FormCommandItem (m, PASTE_MENU_ITEM, NULL, VIB_MSG_PASTE);
  deleteItem = FormCommandItem (m, CLEAR_MENU_ITEM, NULL, VIB_MSG_DELETE);
  SeparatorItem (m);
  if (extraServices) {
    mssgdup = RegisterFormMenuItemName ("SequinDuplicateItem");
    duplicateItem = FormCommandItem (m, "Duplicate...", NULL, mssgdup);
    SeparatorItem (m);
  }
  if (genomeCenter != NULL || indexerVersion) {
    SetupEditSecondary (m, NULL);
    SeparatorItem (m);
  }
  mssgseq = RegisterFormMenuItemName ("SequinEditSequenceItem");
  mssgalign = RegisterFormMenuItemName ("SequinEditAlignmentItem");
  mssgsub = RegisterFormMenuItemName ("SequinEditSubmitterItem");
  mssgupd = RegisterFormMenuItemName ("SequinUpdateSeqSubmenu");
  if (indexerVersion)
  {
    mssgupd_idx = RegisterFormMenuItemName ("SequinUpdateSeqSubmenuIndexer");
  }
  mssgext = RegisterFormMenuItemName ("SequinExtendSeqSubmenu");
  mssgfeatprop = RegisterFormMenuItemName ("SequinFeaturePropagate");
  mssgadd = RegisterFormMenuItemName ("SequinAddSeqSubmenu");
  editsequenceitem = FormCommandItem (m, "Edit Sequence...", NULL, mssgseq);
  editseqalignitem = FormCommandItem (m, "Alignment Assistant...", NULL, mssgalign);
  editseqalignitem = FormCommandItem (m, "Alignment Assistant...", NULL, mssgalign);
  editseqsubitem = FormCommandItem (m, "Edit Submitter Info...", NULL, mssgsub);
  if (indexerVersion) {
    SeparatorItem (m);
    edithistoryitem = CommandItem (m, "Edit History....", EditSequenceHistory);
  }
  SeparatorItem (m);
  
  if (indexerVersion) {
    updateSeqMenuIndexer = SubMenu (m, "Indexer Update Sequence");
    SetFormMenuItem (NULL, mssgupd_idx, (IteM) updateSeqMenuIndexer);
    CommandItem (updateSeqMenuIndexer, "Single Sequence", TestUpdateSequenceIndexer);
    if (useEntrez) {
      CommandItem (updateSeqMenuIndexer, "Download Accession...", UpdateSequenceViaDownloadIndexer);
    }
    SeparatorItem (updateSeqMenuIndexer);
    CommandItem (updateSeqMenuIndexer, "Multiple Sequences", TestUpdateSequenceSetIndexer);
    
    updateSeqMenu = SubMenu (m, "Public Update Sequence");
    SetFormMenuItem (NULL, mssgupd, (IteM) updateSeqMenu);
    CommandItem (updateSeqMenu, "Single Sequence", TestUpdateSequenceSubmitter);
    if (useEntrez) {
      CommandItem (updateSeqMenu, "Download Accession...", UpdateSequenceViaDownloadSubmitter);
    }
    SeparatorItem (updateSeqMenu);
    CommandItem (updateSeqMenu, "Multiple Sequences", TestUpdateSequenceSetSubmitter);
  }
  else
  {
    updateSeqMenu = SubMenu (m, "Update Sequence");
    SetFormMenuItem (NULL, mssgupd, (IteM) updateSeqMenu);
    CommandItem (updateSeqMenu, "Single Sequence", TestUpdateSequenceSubmitter);
    if (useEntrez) {
      CommandItem (updateSeqMenu, "Download Accession...", UpdateSequenceViaDownloadSubmitter);
    }
    SeparatorItem (updateSeqMenu);
    CommandItem (updateSeqMenu, "Multiple Sequences", TestUpdateSequenceSetSubmitter);
  }
  
  extendSeqMenu = SubMenu (m, "Extend Sequence");
  SetFormMenuItem (NULL, mssgext, (IteM) extendSeqMenu);
  CommandItem (extendSeqMenu, "Read FASTA File...", ExtendSeqWithFASTA);
  CommandItem (extendSeqMenu, "Read Sequence Record...", ExtendSeqWithRec);
  CommandItem (extendSeqMenu, "Read FASTA or ASN.1 Set...", ExtendFastaSet);
  if (useEntrez) {
    CommandItem (extendSeqMenu, "Download Accession...", ExtendSeqWithAcc);
  }
  SeparatorItem (m);
  featPropItem = CommandItem (m, "Feature Propagate...", NewFeaturePropagate);
  SetFormMenuItem (NULL, mssgfeatprop, featPropItem);
  SeparatorItem (m);
  addSeqMenu = SubMenu (m, "Add Sequence");
  SetFormMenuItem (NULL, mssgadd, (IteM) addSeqMenu);
  CommandItem (addSeqMenu, "Add FASTA File...", AddSeqWithFASTA);
  CommandItem (addSeqMenu, "Add Sequence Record...", AddSeqWithRec);
  if (! extraServices) {
    SeparatorItem (m);
    parseFileItem = CommandItem (m, "Parse File to Source", ParseFileToSource);
  }

  m = PulldownMenu (NULL, "Search");
  findItem = CommandItem (m, "Find ASN.1.../F", FindStringProc);
  findFFItem = CommandItem (m, "Find FlatFile.../G", FindFlatfileProc);
  SeparatorItem (m);
  findGeneItem = CommandItem (m, "Find by Gene...", FindGeneProc);
  findProtItem = CommandItem (m, "Find by Protein...", FindProtProc);
  findPosItem = CommandItem (m, "Find by Position...", FindPosProc);
  SeparatorItem (m);
  if (indexerVersion) {
    validateMenu = SubMenu (m, "Validate");
    CommandItem (validateMenu, "Validate Record/ V", ValSeqEntryProc);
    CommandItem (validateMenu, "Validate no Alignments", ValSeqEntryProcNoAln);
    CommandItem (validateMenu, "Validate check Inference", ValSeqEntryProcInfAccn);
    SeparatorItem (validateMenu);
    CommandItem (validateMenu, "Validate Inst", ValSeqEntryProcInst);
    CommandItem (validateMenu, "Validate Hist", ValSeqEntryProcHist);
    CommandItem (validateMenu, "Validate Context", ValSeqEntryProcContext);
    CommandItem (validateMenu, "Validate Graph", ValSeqEntryProcGraph);
    CommandItem (validateMenu, "Validate Set", ValSeqEntryProcSet);
    CommandItem (validateMenu, "Validate Feat", ValSeqEntryProcFeat);
    CommandItem (validateMenu, "Validate Desc", ValSeqEntryProcDesc);
  } else {
    validateItem = CommandItem (m, "Validate", ValSeqEntryProc);
  }
#ifdef USE_SPELL
  SeparatorItem (m);
  spellItem = CommandItem (m, "Spell Check...", SpellCheckSeqEntryProc);
#endif
/*#ifdef USE_BLAST*/
  if (useBlast) {
    SeparatorItem (m);
    if (indexerVersion) {
      cddBlastItem = CommandItem (m, "CDD BLAST...", SimpleCDDBlastProc);
    }
    if (extraServices) {
      cddSearchMenu = SubMenu (m, "CDD Search");
      CommandItem (cddSearchMenu, "Features", SimpleCDDSearchFeatProc);
      CommandItem (cddSearchMenu, "Alignments", SimpleCDDSearchAlignProc);
    } else {
      cddSearchItem = CommandItem (m, "CDD Search", SimpleCDDSearchFeatProc);
    }
  }
/*#endif*/
  SeparatorItem (m);
  vecscreenMenu = SubMenu (m, "Vector Screen");
  CommandItem (vecscreenMenu, "UniVec", SimpleUniVecScreenProc);
  CommandItem (vecscreenMenu, "UniVec Core", SimpleUniVecCoreScreenProc);
  SeparatorItem (m);
  orfItem = CommandItem (m, "ORF Finder...", FindOrf);
  /*
  aluItem = CommandItem (m, "Repeat Finder...", FindAlu);
  */
  SeparatorItem (m);
  targetItem = CommandItem (m, "Select Target...", DoChangeTarget);

  m = PulldownMenu (NULL, "Options");
  if (useEntrez) {
    preferencesItem = CommandItem (m, "Preferences...", PreferencesProc);
    SeparatorItem (m);
  }
  /*
  sub = SubMenu (m, "Font Selection");
  if (useEntrez) {
    docsumfontItem = CommandItem (sub, "DocSum Font...", DocSumFontChangeProc);
  }
  displayfontItem = CommandItem (sub, "Display Font...", DisplayFontChangeProc);
  SeparatorItem (m);
  legendItem = CreateLegendItem (m, NULL);
  SeparatorItem (m);
  sub = SubMenu (m, "Query Style");
  queryChoice = CreateQueryTypeChoice (sub, NULL);
  clearUnusedItem = CreateClearUnusedItem (m, NULL);
  SeparatorItem (m);
  sub = SubMenu (m, "Neighbor Policy");
  neighborChoice = CreateNeighborDelayChoice (sub, NULL);
  SeparatorItem (m);
  LoadDocsumOptionsMenu (m);
  seqviewprocs.alignWithChecked = entrezglobals.alignWithChecked;
  seqviewprocs.alignDefault = entrezglobals.alignDefault;
  */
  if (indexerVersion) {
    /*
    SeparatorItem (m);
    */
    oldAsnItem = StatusItem (m, "Use Old ASN.1", ToggleOldAsnProc);
    SetStatus (oldAsnItem, leaveAsOldAsn);
  }

/*#ifdef EXTRA_SERVICES*/
  if (extraServices) {
    specialMenu = PulldownMenu (NULL, "Special");
    SetupSpecialMenu (specialMenu, NULL);
    projectsMenu = PulldownMenu (NULL, "Projects");
    MakeSpecialProjectsMenu (projectsMenu, NULL);

  }
/*#endif*/

  m = PulldownMenu (NULL, "Misc");
  /*
  CommandItem (m, "Style Manager...", StyleManagerProc);
  SeparatorItem (m);
  */
  CommandItem (m, "Net Configure...", NetConfigureProc);
  if (useEntrez) {
    /*
    SeparatorItem (m);
    CommandItem (m, "Entrez2 Query...", Entrez2QueryProc);
    */
/*
#ifndef WIN16
    if (BiostrucAvail ()) {
      SeparatorItem (m);
      CommandItem (m, "Cn3D Window...", Cn3DWinShowProc);
    }
#endif
*/
  }
  if (useDesktop) {
    SeparatorItem (m);
    VSMAddToMenu (m, VSM_DESKTOP);
  }

  /*
  if (indexerVersion) {
    SeparatorItem (m);
    CommandItem (m, "Object Manager Report", ObjMgrReport);
  }
  */

  analysisMenu = CreateAnalysisMenu (NULL, NULL, TRUE, TRUE);

  newFeatMenu = PulldownMenu (NULL, "Annotate");
  SetupNewFeaturesMenu (newFeatMenu, NULL);
  SeparatorItem (newFeatMenu);
  batchApplyMenu = SubMenu (newFeatMenu, "Batch Feature Apply");
  SetupBatchApplyMenu (batchApplyMenu, NULL);
  batchEditMenu = SubMenu (newFeatMenu, "Batch Feature Edit");
  SetupBatchEditMenu (batchEditMenu, NULL);
  CommandItem (newFeatMenu, "Batch Apply Molecule Type", ExternalApplyMoleculeType);
  CommandItem (newFeatMenu, "Set Release Date", SetReleaseDate);
  SeparatorItem (newFeatMenu);
  CommandItem (newFeatMenu, "ORF Finder", FindOrf);
  CommandItem (newFeatMenu, "Import Source Table", ParseFileToSource);
  SeparatorItem (newFeatMenu);
  newPubMenu = SubMenu (newFeatMenu, "Publications");
  SetupNewPublicationsMenu (newPubMenu, NULL);
  SeparatorItem (newFeatMenu);
  newDescMenu = SubMenu (newFeatMenu, "Descriptors");
  SetupNewDescriptorsMenu (newDescMenu, NULL);
  SeparatorItem (newFeatMenu);
  advTableMenu = SubMenu (newFeatMenu, "Advanced Table Readers");
  CommandItem (advTableMenu, "Load Structured Comments from Table", SubmitterCreateStructuredComments);
  sucItem = CommandItem (newFeatMenu, "Sort Unique Count By Group", SUCSubmitterProc);
}
#endif

extern Boolean allowUnableToProcessMessage;
Boolean allowUnableToProcessMessage = TRUE;

static CharPtr canfeatpropagate =
"Since this record has an alignment, you can annotate features on " \
"one record and then use use feature propagation to annotate the "\
"equivalent features on the other records. This is much faster " \
"than annotating everything manually.";

extern Int2 LIBCALLBACK AssemblyUserGenFunc (Pointer data);
static void s_GetTpaInfo (SequencesFormPtr sqfp)
{

  ObjMgrPtr      omp;
  OMProcControl  ompc;
  ObjMgrProcPtr  ompp;
  Int2           retval;
  
  omp = ObjMgrGet ();
  ompp = ObjMgrProcFind (omp, 0, "Edit Assembly User Desc", 0);
  MemSet ((Pointer) &ompc, 0, sizeof (OMProcControl));
  ompc.input_entityID = 1;
  ompc.input_itemID = 1;
  ompc.input_itemtype = OBJ_BIOSEQ;
  ompc.proc = ompp;
  
  retval = AssemblyUserGenFunc (&ompc);
  if (retval == OM_MSG_RET_ERROR) {
    ErrShow ();
  }
  Update ();
}

static CharPtr  tpaString = NULL;

static void FinishPuttingTogether (ForM f)

{
  BaseFormPtr       bfp;
  BioseqSetPtr      bssp;
  Uint2             entityID = 0;
  Int2              handled;
  ObjMgrDataPtr     omdp;
  SubmitBlockPtr    sbp;
  SeqEntryPtr       sep = NULL;
  SequencesFormPtr  sqfp;
  SeqSubmitPtr      ssp;

  bfp = (BaseFormPtr) GetObjectExtra (f);
  if (bfp != NULL) {
    sep = (SeqEntryPtr) FormToPointer (bfp->form);
    if (sep != NULL) {
/*#ifdef USE_TAXON*/
      if (! leaveAsOldAsn) {
        MySeqEntryToAsn3 (sep, TRUE, FALSE, FALSE);
      }
/*#endif*/
      entityID = PackageFormResults (globalsbp, sep, TRUE);
      sqfp = (SequencesFormPtr) bfp;
      if (SEQ_TPA_SUBMISSION == sqfp->submType && entityID > 0) {
        omdp = ObjMgrGetData (entityID);
        if (omdp != NULL && omdp->datatype == OBJ_SEQSUB) {
          ssp = (SeqSubmitPtr) omdp->dataptr;
          if (ssp != NULL && ssp->datatype == 1) {
            sbp = ssp->sub;
            if (sbp != NULL) {
              if (sbp->comment == NULL && StringDoesHaveText (tpaString)) {
                sbp->comment = tpaString;
                tpaString = NULL;
              }
            }
          }
        }
      }
      globalsbp = NULL;
      WatchCursor ();
      seqviewprocs.forceSeparateViewer = TRUE;
      SeqEntrySetScope (NULL);
      handled = GatherProcLaunch (OMPROC_VIEW, FALSE, entityID, 1,
                                  OBJ_BIOSEQ, 0, 0, OBJ_BIOSEQ, 0);
      ArrowCursor ();
      if (handled != OM_MSG_RET_DONE || handled == OM_MSG_RET_NOPROC) {
        Message (MSG_FATAL, "Unable to launch viewer.");
      } else {
        SendHelpScrollMessage (helpForm, "Editing the Record", NULL);
      }
      ObjMgrSetOptions (OM_OPT_FREE_IF_NO_VIEW, entityID);
      ObjMgrSetDirtyFlag (entityID, TRUE);
    } else if (allowUnableToProcessMessage) {
      Message (MSG_FATAL, "Unable to process Seq-entry.");
    }
    sqfp = (SequencesFormPtr) bfp;
    if (SEQ_TPA_SUBMISSION == sqfp->submType) {
      s_GetTpaInfo (sqfp);
    }
    /*SetChecklistValue (checklistForm, 5);*/
    Remove (bfp->form);
    if (sep != NULL && entityID > 0 && IS_Bioseq_set (sep)) {
      bssp = (BioseqSetPtr) sep->data.ptrvalue;
      if (bssp != NULL && (IsPopPhyEtcSet (bssp->_class))) {
        if (SeqEntryHasAligns (entityID, sep)) {
          Message (MSG_OK, "%s", canfeatpropagate);
        }
      }
    }
  }
}

static CharPtr multcomponent = "\
ERROR - You may not enter multiple segments for a single sequence submission.\n\
You should either clear the nucleotide and import a single FASTA record, or \n\
return to the Sequence Format form and choose the proper submission type.";


static Boolean SequencesFormHasNoAnnotation (SequencesFormPtr sqfp)
{
  Int2 annot_type;
  Boolean rval = TRUE;

  if (sqfp == NULL) return FALSE;

  annot_type = GetValue (sqfp->annotType);
  
  switch (annot_type) {
    case 1:
      if (!TextHasNoText (sqfp->geneName)) {
        rval = FALSE;
      }
      break;
    case 2:
      if (!TextHasNoText (sqfp->protOrRnaName)) {
        rval = FALSE;
      }
      break;
    case 3:
      if (!TextHasNoText (sqfp->protOrRnaName) || !TextHasNoText (sqfp->protDesc)) {
        rval = FALSE;
      }
      break;
  }
  return rval;
}


/*---------------------------------------------------------------------*/
/*                                                                     */
/* PutItTogether () --                                                 */
/*                                                                     */
/*---------------------------------------------------------------------*/

static void PutItTogether (ButtoN b)

{
  MsgAnswer    ans;
  BaseFormPtr  bfp;
  ValNodePtr   head;
  SeqEntryPtr  prot_list = NULL, nuc_list = NULL;
  SequencesFormPtr sqfp;

  bfp = (BaseFormPtr) GetObjectExtra (b);
  if (bfp != NULL) {
    head = TestForm (bfp->form);
    if (head != NULL) {
      if (! DisplayTestResults (head)) {
        ValNodeFreeData (head);
        return;
      }
      ValNodeFreeData (head);
    }
    if (SequencesFormHasTooManyNucleotides (bfp->form)) {
      Message (MSG_OK, "%s", multcomponent);
      return;
    }

    if (TRUE == HasZeroLengthSequence (bfp->form)) {
      Message (MSG_POSTERR, "One or more of the submitted sequences are "
	       "zero length.  Zero-length sequences are not allowed. Your "
	       "FASTA sequence may be missing a carriage return after the "
	       "definition line, or may have other formatting problems. "
	       "Please check your FASTA file.");
      return;
    }
    
    /* check for proteins, create nucleotide-protein association */
    prot_list = GetSequencesFormProteinList (bfp->form);
    if (prot_list == NULL && ! SequencesFormHasProteins (bfp->form)) {
      sqfp = (SequencesFormPtr) GetObjectExtra (bfp->form);
      if (SequencesFormHasNoAnnotation (sqfp)) {
        ans = Message (MSG_OKC, "You have not entered proteins and have not created any features.  Is this correct?");
      } else {
        ans = Message (MSG_OKC, "You have not entered proteins.  Is this correct?");
      }
      if (ans == ANS_CANCEL) return;
    }
    else if (prot_list != NULL)
    {
      nuc_list = GetSequencesFormNucleotideList (bfp->form);
      sqfp = (SequencesFormPtr) GetObjectExtra (bfp->form);
      if (sqfp != NULL) 
      {     
        sqfp->nuc_prot_assoc_list = FreeAssociationList (sqfp->nuc_prot_assoc_list);
        sqfp->nuc_prot_assoc_list = AssignProteinsForSequenceSet (nuc_list, prot_list, FALSE);
        if (sqfp->nuc_prot_assoc_list == NULL)
        {
          return;
        }
        
        
        
      }
    }

    Hide (bfp->form);
    ConfirmSequencesFormParsing (bfp->form, FinishPuttingTogether);
  }
}

static void GetOrgAndSeq (ButtoN b);
static void BackToStartup (ButtoN b);
static void BackToFormat (ButtoN b)

{
  MsgAnswer    ans;
  BaseFormPtr  bfp;

  bfp = (BaseFormPtr) GetObjectExtra (b);
  if (bfp != NULL) {
    ans = Message (MSG_OKC, "Are you sure?  Organism and sequence information will be lost.");
    if (ans == ANS_CANCEL) return;
    Remove (bfp->form);
    WatchCursor ();
    Update ();
    PointerToForm (formatForm, &globalFormatBlock);
    Show (formatForm);
    Select (formatForm);
    SendHelpScrollMessage (helpForm, "Sequence Format Form", NULL);
    ArrowCursor ();
    Update ();
  }
}

static void FinishOrgAndSeq (void)

{
  MonitorPtr  mon;
  ForM        w;

  WatchCursor ();
  mon = MonitorStrNewEx ("Sequin New Submission", 30, FALSE);
  MonitorStrValue (mon, "Creating Sequences Form");
  Update ();
  w = CreateInitOrgNucProtForm (-5, -67, "Organism and Sequences",
                                &globalFormatBlock,
                                PutItTogether, BackToFormat,
                                OrgAndSeqsActivateProc);
  ArrowCursor ();
  /*SetChecklistValue (checklistForm, 1);*/
  MonitorFree (mon);
  Update ();
  if (w != NULL) {
    Show (w);
    Select (w);
    if (globalFormatBlock.seqFormat == SEQ_FMT_ALIGNMENT) {
      SendHelpScrollMessage (helpForm, "Nucleotide Page", "Nucleotide Page for Aligned Data Formats");
    } else {
      SendHelpScrollMessage (helpForm, "Nucleotide Page", "Nucleotide Page for FASTA Data Format");
    }
  } else {
    Message (MSG_FATAL, "Unable to create window.");
  }
  Update ();
}

static void BackToSubmitter (ButtoN b)

{
  MsgAnswer    ans;

  ans = Message (MSG_OKC, "Are you sure?  Format information will be lost.");
  if (ans == ANS_CANCEL) return;
  Hide (formatForm);
  Update ();
  PointerToForm (initSubmitForm, globalsbp);
  globalsbp = SequinBlockFree (globalsbp);
  Show (initSubmitForm);
  Select (initSubmitForm);
  SendHelpScrollMessage (helpForm, "Submitting Authors Form", NULL);
  Update ();
  globalFormatBlock.seqPackage = SEQ_PKG_SINGLE;
  globalFormatBlock.seqFormat = SEQ_FMT_FASTA;
  globalFormatBlock.numSeqs = 0;
  globalFormatBlock.submType = SEQ_ORIG_SUBMISSION;
}

static void GetFormat (ButtoN b)

{
  ValNodePtr   head;

  head = TestForm (initSubmitForm);
  if (head != NULL) {
    if (! DisplayTestResults (head)) {
      ValNodeFreeData (head);
      return;
    }
    ValNodeFreeData (head);
  }
  Hide (initSubmitForm);
  globalsbp = (SequinBlockPtr) FormToPointer (initSubmitForm);
  if (globalsbp == NULL) {
    Message (MSG_OK, "Record will be a Seq-entry instead of a Seq-submit.");
  }
  PointerToForm (formatForm, &globalFormatBlock);
  Show (formatForm);
  Select (formatForm);
  SendHelpScrollMessage (helpForm, "Sequence Format Form", NULL);
  Update ();
}

static WindoW   tpaWindow = NULL;
static TexT     tpaText = NULL;
static ButtoN   tpaNext = NULL;
/* tpaString defined above FinishPuttingTogether */

static void TpaPrev (ButtoN b)

{
  Hide (tpaWindow);
  tpaString = MemFree (tpaString);
  SetTitle (tpaText, "");
  Show (formatForm);
  Select (formatForm);
  SendHelpScrollMessage (helpForm, "Sequence Format Form", NULL);
  Update ();
}

static void TpaNext (ButtoN b)

{
  tpaString = MemFree (tpaString);
  tpaString = SaveStringFromText (tpaText);
  if (StringHasNoText (tpaString)) {
    Message (MSG_OK, "The requested information is required in order for you to be able to proceed with a TPA submission");
    return;
  }
  Hide (tpaWindow);
  WatchCursor ();
  FinishOrgAndSeq ();
}

static void TpaText (TexT t)

{
  if (TextHasNoText (t)) {
    SafeDisable (tpaNext);
  } else {
    SafeEnable (tpaNext);
  }
}

static CharPtr tpaMssg = "\
Third party annotation records require a publication describing the biological \
experiments used as evidence for the annotation.  Please provide information \
regarding the nature of these experiments.";

static void DoTpaForm (void)

{
  GrouP   c, h, p;

  if (tpaWindow == NULL) {
    tpaWindow = FixedWindow (-50, -33, -10, -10, "TPA Evidence", NULL);
    h = HiddenGroup (tpaWindow, -1, 0, NULL);
    SetGroupSpacing (h, 10, 10);

    p = MultiLinePrompt (h, tpaMssg, 30 * stdCharWidth, programFont);

    tpaText = ScrollText (h, 30, 5, programFont, TRUE, TpaText);

    c = HiddenGroup (h, 2, 0, NULL);
    PushButton (c, "<< Prev Form", TpaPrev);
    tpaNext = PushButton (c, "Next Form >>", TpaNext);

    AlignObjects (ALIGN_CENTER, (HANDLE) p, (HANDLE) tpaText, (HANDLE) c, NULL);

    RealizeWindow (tpaWindow);
  }
  tpaString = MemFree (tpaString);
  SafeSetTitle (tpaText, "");
  SafeDisable (tpaNext);
  Show (tpaWindow);
  Select (tpaWindow);
}

static void GetOrgAndSeq (ButtoN b)

{
  /*
  MsgAnswer       ans;
  */
  FormatBlockPtr  fbp;
  Boolean         is_tpa = FALSE;

  Hide (formatForm);
  fbp = (FormatBlockPtr) FormToPointer (formatForm);
  if (fbp != NULL) {
    globalFormatBlock.seqPackage = fbp->seqPackage;
    globalFormatBlock.seqFormat = fbp->seqFormat;
    globalFormatBlock.numSeqs = fbp->numSeqs;
    globalFormatBlock.submType = fbp->submType;
    is_tpa = (Boolean) (globalFormatBlock.submType == SEQ_TPA_SUBMISSION);
  }
  MemFree (fbp);
  if (is_tpa) {
    DoTpaForm ();
    /*
    ans = Message (MSG_YN, "%s", tpaMssg);
    if (ans == ANS_YES) {
      WatchCursor ();
      FinishOrgAndSeq ();
    } else {
      Show (formatForm);
      Select (formatForm);
      SendHelpScrollMessage (helpForm, "Sequence Format Form", NULL);
      Update ();
    }
    */
  } else {
    WatchCursor ();
    FinishOrgAndSeq ();
  }
}

static void BackToStartup (ButtoN b)

{
  MsgAnswer    ans;

  ans = Message (MSG_OKC, "Are you sure?  Submitter information will be lost.");
  if (ans == ANS_CANCEL) return;
  Hide (initSubmitForm);
  Update ();
  Show (startupForm);
  Select (startupForm);
  SendHelpScrollMessage (helpForm, "Introduction", NULL);
  Update ();
}

static void ShowHelp (ButtoN b)

{
  if (helpForm == NULL) {
    WatchCursor ();
    helpForm = CreateHelpForm (-95, -5, "Sequin Help", "sequin.hlp",
                               HideHelpForm, HelpActivateProc);
    ArrowCursor ();
  }
  if (helpForm != NULL) {
    Show (helpForm);
    Select (helpForm);
  }
}

static void StartNew (ButtoN b)

{
  Hide (startupForm);
  Update ();
  PointerToForm (initSubmitForm, NULL);
  Show (initSubmitForm);
  Select (initSubmitForm);
  SendHelpScrollMessage (helpForm, "Submitting Authors Form", NULL);
  Update ();
}

static void CloseSubmissionTemplateEditor (Pointer userdata, WindoW w)
{
  Remove (w);
  Show (startupForm);
  Update ();  
}

static void CreateSubmissionTemplate (ButtoN b)

{
  WindoW w;
  
  Hide (startupForm);
  Update ();
  
  w = (WindoW) CreateSubmitTemplateEditorForm (-50, -33, "Submission Template Editor",
                                      CloseSubmissionTemplateEditor, NULL);
  
  Show (w); 
}

static void FinishGenome (ButtoN b)

{
  BaseFormPtr   bfp;
  Uint2         entityID;
  Int2          handled;
  SeqEntryPtr   sep;
  SeqSubmitPtr  ssp;

  bfp = (BaseFormPtr) GetObjectExtra (b);
  if (bfp == NULL) return;
  ssp = (SeqSubmitPtr) FormToPointer (bfp->form);
  if (ssp == NULL) return;
  if (ssp->datatype != 1) return;
  sep = (SeqEntryPtr) ssp->data;
  if (sep == NULL) return;
  ObjMgrConnect (OBJ_SEQENTRY, sep->data.ptrvalue, OBJ_SEQSUB, (Pointer) ssp);
  if (! ObjMgrRegister (OBJ_SEQSUB, (Pointer) ssp)) {
    ErrPostEx (SEV_ERROR, 0, 0, "ObjMgrRegister failed.");
  }
  if (! leaveAsOldAsn) {
    MySeqEntryToAsn3 (sep, TRUE, FALSE, FALSE);
  }
  entityID = ObjMgrGetEntityIDForChoice (sep);
  WatchCursor ();
  seqviewprocs.forceSeparateViewer = TRUE;
  SeqEntrySetScope (NULL);
  handled = GatherProcLaunch (OMPROC_VIEW, FALSE, entityID, 1,
                              OBJ_BIOSEQ, 0, 0, OBJ_BIOSEQ, 0);
  ArrowCursor ();
  if (handled != OM_MSG_RET_DONE || handled == OM_MSG_RET_NOPROC) {
    Message (MSG_FATAL, "Unable to launch viewer.");
  } else {
    SendHelpScrollMessage (helpForm, "Editing the Record", NULL);
  }
  ObjMgrSetOptions (OM_OPT_FREE_IF_NO_VIEW, entityID);
  ObjMgrSetDirtyFlag (entityID, TRUE);
}

static void CancelGenome (ButtoN b)

{
  Show (startupForm);
  Select (startupForm);
  SendHelpScrollMessage (helpForm, "Introduction", NULL);
  Update ();
}

static void StartFa2htgs (ButtoN b)

{
  ForM  f;

  Hide (startupForm);
  Update ();
  f = CreateGenomeCenterForm (-50, -33, "HTGS Submission",
                              FinishGenome, CancelGenome, FALSE, FALSE,
                              GenomeFormActivateProc);
  if (f != NULL) {
    Show (f);
    Select (f);
  } else {
    Show (startupForm);
  }
  Update ();
}

static void StartPhrap (ButtoN b)

{
  ForM  f;

  Hide (startupForm);
  Update ();
  f = CreateGenomeCenterForm (-50, -33, "PHRAP Submission",
                              FinishGenome, CancelGenome, TRUE, FALSE,
                              GenomeFormActivateProc);
  if (f != NULL) {
    Show (f);
    Select (f);
  } else {
    Show (startupForm);
  }
  Update ();
}

static void BuildContig (ButtoN b)

{
  ForM  f;

  Hide (startupForm);
  Update ();
  f = CreateGenomeCenterForm (-50, -33, "Contig Submission",
                              FinishGenome, CancelGenome, FALSE, TRUE,
                              GenomeFormActivateProc);
  if (f != NULL) {
    Show (f);
    Select (f);
  } else {
    Show (startupForm);
  }
  Update ();
/*
  Hide (startupForm);
  Update ();
  if (! DoBuildContig ()) {
    Show (startupForm);
  }
  Update ();
*/
}

static void ReadOld (ButtoN b)

{
  Hide (startupForm);
  Update ();
  if (! CommonReadNewAsnProc (NULL, FALSE, TRUE)) {
    Show (startupForm);
    Select (startupForm);
  }
}

static void ReadSettings (void)

{
  Char  str [64];

  useDesktop = FALSE;
  useEntrez = FALSE;
  useLocal = FALSE;
  useIdLookup = FALSE;
  useBlast = FALSE;
  useMedarch = FALSE;
  newMedarch = FALSE;
  useTaxon = FALSE;
  allowDownload = FALSE;
  extraServices = FALSE;
  indexerVersion = FALSE;

  genomeCenter = NULL;

#ifdef INTERNAL_NCBI_SEQUIN
  indexerVersion = TRUE;
  extraServices = TRUE;
  useDesktop = TRUE;
  useEntrez = TRUE;
  useLocal = TRUE;
  useIdLookup = TRUE;
  useBlast = TRUE;
  useMedarch = TRUE;
  newMedarch = TRUE;
  useTaxon = TRUE;
  allowDownload = TRUE;
#endif
  if (GetSequinAppParam ("SETTINGS", "INDEXERVERSION", NULL, str, sizeof (str))) {
    if (StringICmp (str, "TRUE") == 0) {
      indexerVersion = TRUE;
      extraServices = TRUE;
      useDesktop = TRUE;
      useEntrez = TRUE;
      useLocal = TRUE;
      useIdLookup = TRUE;
      useBlast = TRUE;
      useMedarch = TRUE;
      newMedarch = TRUE;
      useTaxon = TRUE;
      allowDownload = TRUE;
    }
  }

#ifdef PUBLIC_NETWORK_SEQUIN
  useDesktop = TRUE;
  useEntrez = TRUE;
  useLocal = TRUE;
  useIdLookup = TRUE;
  useBlast = TRUE;
  useMedarch = TRUE;
  newMedarch = TRUE;
  useTaxon = TRUE;
  allowDownload = TRUE;
#endif
  if (GetSequinAppParam ("SETTINGS", "PUBLICNETWORKSEQUIN", NULL, str, sizeof (str))) {
    if (StringICmp (str, "TRUE") == 0) {
      useDesktop = TRUE;
      useEntrez = TRUE;
      useLocal = TRUE;
      useIdLookup = TRUE;
      useBlast = TRUE;
      useMedarch = TRUE;
      newMedarch = TRUE;
      useTaxon = TRUE;
      allowDownload = TRUE;
    }
  }

#ifdef USE_DESKTOP
  useDesktop = TRUE;
#endif
  if (GetSequinAppParam ("SETTINGS", "USEDESKTOP", NULL, str, sizeof (str))) {
    if (StringICmp (str, "TRUE") == 0) {
      useDesktop = TRUE;
    }
  }

#ifdef USE_ENTREZ
  useEntrez = TRUE;
#endif
  if (GetSequinAppParam ("SETTINGS", "USEENTREZ", NULL, str, sizeof (str))) {
    if (StringICmp (str, "TRUE") == 0) {
      useEntrez = TRUE;
    }
  }

#ifdef USE_BLAST
  useBlast = TRUE;
#endif
  if (GetSequinAppParam ("SETTINGS", "USEBLAST", NULL, str, sizeof (str))) {
    if (StringICmp (str, "TRUE") == 0) {
      useBlast = TRUE;
    }
  }

#ifdef USE_MEDARCH
  useMedarch = TRUE;
  newMedarch = TRUE;
#endif
  if (GetSequinAppParam ("SETTINGS", "USEMEDARCH", NULL, str, sizeof (str))) {
    if (StringICmp (str, "TRUE") == 0) {
      useMedarch = TRUE;
      newMedarch = TRUE;
    }
  }
  if (GetSequinAppParam ("SETTINGS", "NEWMEDARCH", NULL, str, sizeof (str))) {
    if (StringICmp (str, "FALSE") == 0) {
      newMedarch = FALSE;
    }
  }

#ifdef USE_TAXON
  useTaxon = TRUE;
#endif
  if (GetSequinAppParam ("SETTINGS", "USETAXON", NULL, str, sizeof (str))) {
    if (StringICmp (str, "TRUE") == 0) {
      useTaxon = TRUE;
    }
  }

#ifdef ALLOW_DOWNLOAD
  allowDownload = TRUE;
#endif
  if (GetSequinAppParam ("SETTINGS", "ALLOWDOWNLOAD", NULL, str, sizeof (str))) {
    if (StringICmp (str, "TRUE") == 0) {
      allowDownload = TRUE;
      useEntrez = TRUE;
    }
  }

#ifdef EXTRA_SERVICES
  extraServices = TRUE;
#endif
  if (GetSequinAppParam ("SETTINGS", "USEEXTRASERVICES", NULL, str, sizeof (str))) {
    if (StringICmp (str, "TRUE") == 0) {
      extraServices = TRUE;
    }
  }

  if (GetSequinAppParam ("SETTINGS", "SUPPRESSENTREZ", NULL, str, sizeof (str))) {
    if (StringICmp (str, "TRUE") == 0) {
      allowDownload = FALSE;
      useEntrez = FALSE;
      useMedarch = FALSE;
      useTaxon = FALSE;
    }
  }

  useSeqFetch = useEntrez;

  if (GetSequinAppParam ("SETTINGS", "SUPPRESSSEQFETCH", NULL, str, sizeof (str))) {
    if (StringICmp (str, "TRUE") == 0) {
      useSeqFetch = FALSE;
    }
  }

  if (GetSequinAppParam ("SETTINGS", "SUPPRESSIDLOOKUP", NULL, str, sizeof (str))) {
    if (StringICmp (str, "TRUE") == 0) {
      useIdLookup = FALSE;
    }
  }

#ifdef USE_LOCAL
  useLocal = TRUE;
#endif
  if (GetSequinAppParam ("SETTINGS", "USELOCAL", NULL, str, sizeof (str))) {
    if (StringICmp (str, "TRUE") == 0) {
      useLocal = TRUE;
    } else if (StringICmp (str, "FALSE") == 0) {
      useLocal = FALSE;
    }
  }

  if (GetSequinAppParam ("SETTINGS", "SUPPRESSLOCAL", NULL, str, sizeof (str))) {
    if (StringICmp (str, "TRUE") == 0) {
      useLocal = FALSE;
    }
  }

  if (GetSequinAppParam ("SETTINGS", "GENOMETAG", NULL, str, sizeof (str))) {
    TrimSpacesAroundString (str);
    if (! StringHasNoText (str)) {
      genomeCenter = StringSave (str);
      extraServices = TRUE;
    }
  } else if (GetSequinAppParam ("SETTINGS", "GENOMECENTER", NULL, str, sizeof (str))) {
    TrimSpacesAroundString (str);
    if (! StringHasNoText (str)) {
      genomeCenter = StringSave (str);
      extraServices = TRUE;
    }
  }

  loadSaveUidListOK = useEntrez;
  if (GetSequinAppParam ("PREFERENCES", "LOADSAVEUIDLIST", NULL, str, sizeof (str))) {
    if (StringICmp (str, "FALSE") == 0) {
      loadSaveUidListOK = FALSE;
    }
  }

  if (GetSequinAppParam ("SETTINGS", "USEOLDGRAPHICAL", NULL, str, sizeof (str))) {
    if (StringICmp (str, "TRUE") == 0) {
      useOldGraphicView = TRUE;
    }
  }

  if (GetSequinAppParam ("SETTINGS", "USEOLDALIGNMENT", NULL, str, sizeof (str))) {
    if (StringICmp (str, "TRUE") == 0) {
      useOldAlignmentView = TRUE;
    }
  }

  if (GetSequinAppParam ("SETTINGS", "USEOLDSEQVIEW", NULL, str, sizeof (str))) {
    if (StringICmp (str, "TRUE") == 0) {
      useOldSequenceView = TRUE;
    }
  }

  /*
  if (GetSequinAppParam ("SETTINGS", "USEUDV", NULL, str, sizeof (str))) {
    if (StringICmp (str, "TRUE") == 0) {
      useUdv = TRUE;
    }
  }
  */

  if (GetSequinAppParam ("SETTINGS", "EDITGBSOURCEQUALS", NULL, str, sizeof (str))) {
    if (StringICmp (str, "TRUE") == 0) {
      SetAppProperty ("EditGBSourceQuals", (void *) 1024);
    }
  }
}

static void DoQuit (ButtoN b)

{
  QuitProgram ();
}

static void CleanupSequin (void)

{
  FiniSequinExtras ();

  ExitMuskStyles ();
  FreeOrganismTable ();
  FreeGeneticCodes ();
  GeneticCodeSingletonFini ();
  FreePrintOptions ();

/*#ifdef USE_LOCAL*/
  if (useLocal) {
    LocalSeqFetchDisable ();
  }
/*#endif*/

/*#ifdef USE_ENTREZ*/
  if (useEntrez) {
    /* EntrezBioseqFetchDisable (); */
    /*
    if (EntrezIsInited ()) {
      EntrezFini ();
    }
    */
    /* ID1BioseqFetchDisable (); */
    PubSeqFetchDisable ();
  }
/*#endif*/

  TransTableFreeAll ();

  ECNumberFSAFreeAll ();
}

static Pointer SubtoolModeAsnTextFileRead (CharPtr filename,
                                           Uint2Ptr datatypeptr, Uint2Ptr entityIDptr)

{
  AsnIoPtr       aip;
  BioseqSetPtr   bssp = NULL, bssp2;
  ObjMgrData     omdata;
  ObjMgrDataPtr  omdptop = NULL;
  Uint2          parenttype = 0;
  Pointer        parentptr = NULL;
  Pointer        ptr = NULL;
  SeqEntryPtr    sep = NULL, sep1, sep2;

  if (filename == NULL) return NULL;
  if (datatypeptr != NULL) *datatypeptr = 0;
  if (entityIDptr != NULL) *entityIDptr = 0;

  aip = AsnIoOpen (filename, "r");
  if (aip == NULL) return NULL;
  aip->scan_for_start = TRUE;
  while ((bssp2 = BioseqSetAsnRead (aip, NULL)) != NULL) {
    if (sep == NULL) {
      sep2 = SeqEntryNew ();
      if (sep2 != NULL) {
        sep2->choice = 2;
        sep2->data.ptrvalue = bssp2;
        SeqMgrSeqEntry (SM_BIOSEQSET, (Pointer) bssp2, sep2);
        sep = sep2;
        bssp = bssp2;
        SaveSeqEntryObjMgrData (sep, &omdptop, &omdata);
        GetSeqEntryParent (sep, &parentptr, &parenttype);
      }
    } else {
      for (sep1 = bssp->seq_set; sep1->next != NULL; sep1 = sep1->next) continue;
      sep1->next = bssp2->seq_set;
      bssp2->seq_set = NULL;
      BioseqSetFree (bssp2);
    }
  }
  if (sep != NULL) {
    SeqMgrLinkSeqEntry (sep, parenttype, parentptr);
    RestoreSeqEntryObjMgrData (sep, omdptop, &omdata);
  }
  AsnIoClose (aip);
  ptr = (Pointer) bssp;
  if (ptr == NULL) {
    ErrPostEx (SEV_ERROR,0,0,"Couldn't read [%s], type [OBJ_BIOSEQSET]", filename);
  } else {
    if (datatypeptr != NULL) *datatypeptr = OBJ_BIOSEQSET;
    if (entityIDptr != NULL) *entityIDptr = ObjMgrRegister (OBJ_BIOSEQSET, ptr);
  }
  return ptr;
}

#ifdef USE_SMARTNET
static CharPtr dirsubfetchproc = "DirSubBioseqFetch";

static CharPtr dirsubfetchcmd = NULL;

extern Pointer ReadFromDirSub (CharPtr accn, Uint2Ptr datatype, Uint2Ptr entityID);
extern Pointer ReadFromDirSub (CharPtr accn, Uint2Ptr datatype, Uint2Ptr entityID)

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
  if (! dirsubMode) return NULL;
  if (StringHasNoText (accn)) return NULL;

  if (dirsubfetchcmd == NULL) {
    if (GetAppParam ("SEQUIN", "DIRSUB", "FETCHSCRIPT", NULL, cmmd, sizeof (cmmd))) {
    	dirsubfetchcmd = StringSaveNoNull (cmmd);
    }
  }
  if (dirsubfetchcmd == NULL) return NULL;

  TmpNam (path);

#ifdef OS_UNIX
  sprintf (cmmd, "csh %s %s > %s", dirsubfetchcmd, accn, path);
  system (cmmd);
#endif
#ifdef OS_MSWIN
  sprintf (cmmd, "%s %s -o %s", dirsubfetchcmd, accn, path);
  RunSilent (cmmd);
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


static Int2 LIBCALLBACK DirSubBioseqFetchFunc (Pointer data)

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

  if (dirsubfetchcmd == NULL) {
    if (GetAppParam ("SEQUIN", "DIRSUB", "FETCHSCRIPT", NULL, cmmd, sizeof (cmmd))) {
    	dirsubfetchcmd = StringSaveNoNull (cmmd);
    }
  }
  if (dirsubfetchcmd == NULL) return OM_MSG_RET_ERROR;

  TmpNam (path);

#ifdef OS_UNIX
  sprintf (cmmd, "csh %s %s > %s", dirsubfetchcmd, tsip->accession, path);
  system (cmmd);
#endif
#ifdef OS_MSWIN
  sprintf (cmmd, "%s %s -o %s", dirsubfetchcmd, tsip->accession, path);
  RunSilent (cmmd);
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

static Boolean DirSubFetchEnable (void)

{
  ObjMgrProcLoad (OMPROC_FETCH, dirsubfetchproc, dirsubfetchproc,
                  OBJ_SEQID, 0, OBJ_BIOSEQ, 0, NULL,
                  DirSubBioseqFetchFunc, PROC_PRIORITY_DEFAULT);
  return TRUE;
}

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
  if (! dirsubMode) return NULL;
  if (StringHasNoText (accn)) return NULL;

  if (smartfetchcmd == NULL) {
    if (GetAppParam ("SEQUIN", "SMART", "FETCHSCRIPT", NULL, cmmd, sizeof (cmmd))) {
    	smartfetchcmd = StringSaveNoNull (cmmd);
    }
  }
  if (smartfetchcmd == NULL) return NULL;

  TmpNam (path);

#ifdef OS_UNIX
  sprintf (cmmd, "csh %s %s > %s", smartfetchcmd, accn, path);
  system (cmmd);
#endif
#ifdef OS_MSWIN
  sprintf (cmmd, "%s %s -o %s", smartfetchcmd, accn, path);
  RunSilent (cmmd);
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
  OMUserDataPtr     omudp;

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
  sprintf (cmmd, "csh %s %s > %s", smartfetchcmd, tsip->accession, path);
  system (cmmd);
#endif
#ifdef OS_MSWIN
  sprintf (cmmd, "%s %s -o %s", smartfetchcmd, tsip->accession, path);
  RunSilent (cmmd);
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

  omudp = ObjMgrAddUserData(ompcp->output_entityID, ompp->procid, OMPROC_FETCH, 0);
  

  return OM_MSG_RET_DONE;
}

static Boolean SmartFetchEnable (void)

{
  ObjMgrProcLoad (OMPROC_FETCH, smartfetchproc, smartfetchproc,
                  OBJ_SEQID, 0, OBJ_BIOSEQ, 0, NULL,
                  SmartBioseqFetchFunc, PROC_PRIORITY_DEFAULT);
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
  if (! dirsubMode) return NULL;
  if (StringHasNoText (accn)) return NULL;

  if (hupfetchcmd == NULL) {
    if (GetAppParam ("SEQUIN", "HUP", "FETCHSCRIPT", NULL, cmmd, sizeof (cmmd))) {
    	hupfetchcmd = StringSaveNoNull (cmmd);
    }
  }
  if (hupfetchcmd == NULL) return NULL;

  TmpNam (path);

#ifdef OS_UNIX
  sprintf (cmmd, "csh %s %s > %s", hupfetchcmd, accn, path);
  system (cmmd);
#endif
#ifdef OS_MSWIN
  sprintf (cmmd, "%s %s -o %s", hupfetchcmd, accn, path);
  RunSilent (cmmd);
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
  OMUserDataPtr     omudp;

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
  sprintf (cmmd, "csh %s %s > %s 2>&1", hupfetchcmd, tsip->accession, path);
  system (cmmd);
#endif
#ifdef OS_MSWIN
  sprintf (cmmd, "%s %s -o %s", hupfetchcmd, tsip->accession, path);
  RunSilent (cmmd);
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

  omudp = ObjMgrAddUserData(ompcp->output_entityID, ompp->procid, OMPROC_FETCH, 0);
  

  return OM_MSG_RET_DONE;
}

static Boolean HUPFetchEnable (void)

{
  ObjMgrProcLoad (OMPROC_FETCH, hupfetchproc, hupfetchproc,
                  OBJ_SEQID, 0, OBJ_BIOSEQ, 0, NULL,
                  HUPBioseqFetchFunc, PROC_PRIORITY_DEFAULT);
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
  if (! dirsubMode) return NULL;
  if (StringHasNoText (accn)) return NULL;

  if (tpasmartfetchcmd == NULL) {
    if (GetAppParam ("SEQUIN", "TPASMART", "FETCHSCRIPT", NULL, cmmd, sizeof (cmmd))) {
    	tpasmartfetchcmd = StringSaveNoNull (cmmd);
    }
  }
  if (tpasmartfetchcmd == NULL) return NULL;

  TmpNam (path);

#ifdef OS_UNIX
  sprintf (cmmd, "csh %s %s > %s", tpasmartfetchcmd, accn, path);
  system (cmmd);
#endif
#ifdef OS_MSWIN
  sprintf (cmmd, "%s %s -o %s", tpasmartfetchcmd, accn, path);
  RunSilent (cmmd);
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
  sprintf (cmmd, "csh %s %s > %s", tpasmartfetchcmd, tsip->accession, path);
  system (cmmd);
#endif
#ifdef OS_MSWIN
  sprintf (cmmd, "%s %s -o %s", tpasmartfetchcmd, tsip->accession, path);
  RunSilent (cmmd);
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

static VoidPtr LIBCALLBACK SmartUserDataFree (VoidPtr Pointer)
{
    SMUserDataPtr sm_user_data = (SMUserDataPtr) Pointer;

    if(sm_user_data != NULL) {
        MemFree(sm_user_data->header);
        MemFree(sm_user_data);
    }

    return NULL;
}

static VoidPtr LIBCALLBACK DumbUserDataFree (VoidPtr Pointer)

{
  return ValNodeFreeData ((ValNodePtr) Pointer);
}

static void SMCancelAllEdit(void)
{
    ObjMgrPtr      omp;
    OMUserDataPtr  omudp;
    Int2           num;
    SMUserDataPtr  sm_usr_data = NULL;
    Uint2 j;

    omp = ObjMgrGet ();
    num = omp->HighestEntityID;
    
    for(j = 1; j <= omp->HighestEntityID; j++) {
        if((omudp = ObjMgrGetUserData(j, 0,0, SMART_KEY)) != NULL) {
            if((sm_usr_data = 
                (SMUserDataPtr) omudp->userdata.ptrvalue) != NULL &&
               sm_usr_data->fd != 0) {
                sm_usr_data->header->status = SMStatClosed;
                SMSendMsgToClient(sm_usr_data);
                SmartUserDataFree(omudp->userdata.ptrvalue);                
                omudp->userdata.ptrvalue = NULL;

                /* Deleting all */
                ObjMgrSendMsg (OM_MSG_DEL, j, 0, 0);
            }
        }
        if ((omudp = ObjMgrGetUserData(j, 0,0, DUMB_KEY)) != NULL) {
          DumbUserDataFree (omudp->userdata.ptrvalue);                
          omudp->userdata.ptrvalue = NULL;
        }
    }
    
    return;
}

extern ForM smartBioseqViewForm;


static void CountBioSources (BioSourcePtr biop, Pointer userdata)
{
  Int4Ptr count;

  if (biop == NULL || userdata == NULL) return;
  count = (Int4Ptr) userdata;
  (*count)++;
}

static Int4 SMReadBioseqObj(VoidPtr data, CharPtr buffer, Int4 length, void* fd)
{
    AsnIoMemPtr    aimp;
    BaseFormPtr    bfp;
    VoidPtr        bio_data;
    Int2           handled; 
    Uint2          entityID = 0;
    OMUserDataPtr  omudp;
    SMMsgHeaderPtr header;
    Int4           headlen;
    BioseqPtr      bsp;
    BioseqSetPtr   bssp = NULL, bssp2;
    SeqEntryPtr    sep = NULL, sep1, sep2;
    ObjMgrData     omdata;
    ObjMgrDataPtr  omdptop = NULL;
    Uint2          parenttype = 0;
    Pointer        parentptr = NULL;
    SMUserDataPtr  sm_user_data;
    Int4           bio_type;
    Int4           num_srcs = 0;
    Boolean        do_taxlookup = TRUE;

    if(buffer == NULL || length < sizeof(SMMsgHeader))
        return -1; 
    
    /* Reading request header */
    
    headlen = sizeof(SMMsgHeader); 
    header = (SMMsgHeaderPtr) MemNew(headlen);
    MemCpy((Pointer) header, buffer, headlen);
    
    /* ----------------------- */

    switch(header->status) {
    case SMTaskEditBinary:
        aimp = AsnIoMemOpen("rb", 
                            (UcharPtr) buffer+headlen, length - headlen);
        break;
    case SMTaskEditText:
        aimp = AsnIoMemOpen("r", (UcharPtr) buffer+headlen, length - headlen);
        break;
    case SMTaskCancel:
        SMCancelAllEdit();
        return TRUE;
    default:
        Message(MSG_ERROR, "This request type (%d) is not implemented yet",
                header->status);
        return FALSE;    
    }

    AsnIoSetErrorMsg(aimp->aip, SMAsnErrorFunc);
    bio_type = header->format;

    switch (header->format) {
    case OBJ_SEQSUB:
        SeqMgrHoldIndexing (TRUE);
        bio_data = (VoidPtr) SeqSubmitAsnRead(aimp->aip, NULL);
        SeqMgrHoldIndexing (FALSE);
        break;
    case OBJ_SEQENTRY:
        SeqMgrHoldIndexing (TRUE);
        bio_data = (VoidPtr) SeqEntryAsnRead(aimp->aip, NULL);
        SeqMgrHoldIndexing (FALSE);
        if((sep = (SeqEntryPtr) bio_data) != NULL) {
            if (sep->choice == 1) {
                bio_type = OBJ_BIOSEQ;
            } else {
                bio_type = OBJ_BIOSEQSET;
            }
        }
        break;
    case OBJ_BIOSEQ:
        SeqMgrHoldIndexing (TRUE);
        bio_data = (VoidPtr) BioseqAsnRead(aimp->aip, NULL);
        SeqMgrHoldIndexing (FALSE);
        bsp = (BioseqPtr) bio_data;
        sep = SeqEntryNew ();
        sep->choice = 1;
        sep->data.ptrvalue = bsp;
        SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) bsp, sep);

        break;
    case OBJ_BIOSEQSET:

        aimp->aip->scan_for_start = TRUE;
        SeqMgrHoldIndexing (TRUE);
        while ((bssp2 = BioseqSetAsnRead (aimp->aip, NULL)) != NULL) {
            if (sep == NULL) {
                sep2 = SeqEntryNew ();
                if (sep2 != NULL) {
                    sep2->choice = 2;
                    sep2->data.ptrvalue = bssp2;
                    SeqMgrSeqEntry (SM_BIOSEQSET, (Pointer) bssp2, sep2);
                    sep = sep2;
                    bssp = bssp2;
                    SaveSeqEntryObjMgrData (sep, &omdptop, &omdata);
                    GetSeqEntryParent (sep, &parentptr, &parenttype);
                }
            } else {
                for (sep1 = bssp->seq_set; sep1->next != NULL; 
                     sep1 = sep1->next) continue;
                sep1->next = bssp2->seq_set;
                bssp2->seq_set = NULL;
                BioseqSetFree (bssp2);
            }
        }
        SeqMgrHoldIndexing (FALSE);
        if (sep != NULL) {
            SeqMgrLinkSeqEntry (sep, parenttype, parentptr);
            RestoreSeqEntryObjMgrData (sep, omdptop, &omdata);
        }
        bio_data = (Pointer) bssp;
        break;
    default:
        Message(MSG_OK, "This datatype: %d is not implemented yet",
                header->format);
        AsnIoMemClose(aimp);  
        return FALSE;
    }
    
    sm_user_data = MemNew(sizeof(SMUserData));
    sm_user_data->fd = (void*) fd;
    sm_user_data->header = header;

    if(bio_data == NULL) {
        Message(MSG_ERROR, "Error parcing or processing BLOB");
        sm_user_data->header->status = SMStatClosed;
        SMSendMsgToClient(sm_user_data); 
        MemFree(sm_user_data->header);
        MemFree(sm_user_data);
    } else {

        if (smartBioseqViewForm != NULL) {
            bfp = (BaseFormPtr) GetObjectExtra (smartBioseqViewForm);
            if (bfp != NULL) {
                /*
                entityID = bfp->input_entityID;
                RemoveSeqEntryViewerEx (bfp->form, FALSE);
                ObjMgrSendMsg (OM_MSG_DEL, entityID, 0, 0);
                bfp->input_entityID = 0;
                */
                seqviewprocs.keepSmartViewerVisible = TRUE;
                SmartnetDoneFunc (bfp);
                seqviewprocs.keepSmartViewerVisible = FALSE;
            }
        }

        entityID = ObjMgrRegister (bio_type, bio_data);
        subtoolEntityID = entityID;

        sep = GetTopSeqEntryForEntityID (entityID);
        VisitBioSourcesInSep (sep, &num_srcs, CountBioSources);
        if (num_srcs > 10000) {
          if (ANS_CANCEL == Message (MSG_OKC, "Record contains %d BioSources.  Do Taxlookup now?", num_srcs)) {
            do_taxlookup = FALSE;
          }
        }
        MySeqEntryToAsn3 (sep, TRUE, FALSE, do_taxlookup);

        /* now instantiating protein titles */
        InstantiateProteinTitles (entityID, NULL);

        if((omudp = ObjMgrGetUserData (entityID, 0, 0, SMART_KEY)) != NULL) {
            SmartUserDataFree(omudp->userdata.ptrvalue);
        } else {
            omudp = ObjMgrAddUserData (entityID, 0, 0, SMART_KEY);
        }

        ObjMgrSetOptions (OM_OPT_FREE_IF_NO_VIEW, entityID);
        
        /* We will set dirty flag from header info */

        if(header->dirty & 0x01) {
            ObjMgrSetDirtyFlag (entityID, TRUE);
        } else {
            ObjMgrSetDirtyFlag (entityID, FALSE);
        }

        if (omudp != NULL) {
            omudp->userdata.ptrvalue = (VoidPtr) sm_user_data;
            omudp->messagefunc = SubtoolModeMsgFunc; /* ? */
            omudp->freefunc = SmartUserDataFree;
        }
 
        if((omudp = ObjMgrGetUserData (entityID, 0, 0, DUMB_KEY)) != NULL) {
            DumbUserDataFree (omudp->userdata.ptrvalue);
        } else {
            omudp = ObjMgrAddUserData (entityID, 0, 0, DUMB_KEY);
        }

        if (omudp != NULL) {
            omudp->userdata.ptrvalue = (VoidPtr) SmartStructureReport (sep);
            omudp->messagefunc = NULL;
            omudp->freefunc = DumbUserDataFree;
        }
 
        seqviewprocs.forceSeparateViewer = FALSE;
        SeqEntrySetScope (NULL);
        handled = GatherProcLaunch (OMPROC_VIEW, FALSE, entityID, 1,
                                    OBJ_BIOSEQ, 0, 0, OBJ_BIOSEQ, 0);
        ArrowCursor ();
        
        if (handled != OM_MSG_RET_DONE || handled == OM_MSG_RET_NOPROC) {
            Message (MSG_FATAL, "Unable to launch viewer.");
            CleanupSequin ();
            return FALSE;
        } 
    }

    AsnIoMemClose(aimp);  
    return TRUE;
}

static Int4 SMWriteBioseqObj(VoidPtr bio_data, 
                             SMUserDataPtr sm_usr_data, 
                             VoidPtr data)
{
    ByteStorePtr    bsp;
    AsnIoBSPtr      aibp;
    CharPtr buffer;
    Int4 length, totlen;

    bsp = BSNew(1024);

    if(sm_usr_data->header->status == SMTaskEditBinary) {  
        aibp = AsnIoBSOpen("wb", bsp);
    } else {
        aibp = AsnIoBSOpen("w", bsp);
    }

    switch (sm_usr_data->header->format) {
    case OBJ_SEQSUB:
        SeqSubmitAsnWrite((SeqSubmitPtr) bio_data, aibp->aip, NULL);
        break;
    case OBJ_SEQENTRY:
        SeqEntryAsnWrite((SeqEntryPtr) bio_data, aibp->aip, NULL);
        break;
    case OBJ_BIOSEQ:
        BioseqAsnWrite((BioseqPtr) bio_data, aibp->aip, NULL);
        break;
    case OBJ_BIOSEQSET:
        BioseqSetAsnWrite((BioseqSetPtr) bio_data, aibp->aip, NULL);
        break;
    default:
        Message(MSG_OK, "This datatype: %d is not implemented yet");
        AsnIoBSClose(aibp);
        BSFree(bsp);
        return FALSE;
    }

    AsnIoBSClose(aibp);
    
    BSSeek(bsp, 0, SEEK_SET);
    
    length = BSLen(bsp);
    buffer = MemNew(length+1);
    BSRead(bsp, buffer, length);
    
    totlen = sizeof(SMMsgHeader) + length;

    if(SMWriteToClient(sm_usr_data->fd, 
                       (CharPtr) &totlen, sizeof(Int4)) != SMNoError) {
        Message(MSG_OK, "Write error. Errno = %d", errno);
        return FALSE;
    }
    
    if(SMWriteToClient(sm_usr_data->fd, (CharPtr) sm_usr_data->header, 
                       sizeof(SMMsgHeader)) != SMNoError) {
        Message(MSG_OK, "Write error. Errno = %d", errno);
        return FALSE;
    }

    if(SMWriteToClient(sm_usr_data->fd, buffer, length) != SMNoError) {
        Message(MSG_OK, "Write error. Errno = %d", errno);
        return FALSE;
    }

    BSFree(bsp);
    MemFree(buffer);
    return TRUE;
}
#endif

static CharPtr comp_months [] = {
  "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", NULL
};

static Boolean MoreThanYearOld (void)

{
  Int4      compday, compmonth, compyear, currday, currmonth, curryear;
  DayTime   dt;
  Int2      i;
  CharPtr   ptr, str;
  Char      tmp [80];
  long int  val;

  if (! GetDayTime (&dt)) return FALSE;
  currmonth = dt.tm_mon + 1;
  currday = dt.tm_mday;
  curryear = dt.tm_year + 1900;

  compmonth = 0;
  compday = 0;
  compyear = 0;

  sprintf (tmp, "%s", date_of_compilation);

  ptr = StringChr (tmp, ' ');
  if (ptr != NULL) {
    *ptr = '\0';
    ptr++;
    if (*ptr == ' ') {
      ptr++;
    }
    for (i = 0; i < 12; i++) {
      if (StringCmp (tmp, comp_months [i]) == 0) {
        compmonth = (Int4) i + 1;
        break;
      }
    }
    str = StringChr (ptr, ' ');
    if (str != NULL) {
      *str = '\0';
      str++;
      if (sscanf (ptr, "%ld", &val) == 1) {
        compday = (Int4) val;
      }
      if (sscanf (str, "%ld", &val) == 1) {
        compyear = (Int4) val;
      }
    }
  }

  if (compmonth == 0 || compyear == 0) return FALSE;

  if (curryear > compyear + 1) return TRUE;
  if (curryear == compyear + 1) {
    if (currmonth > compmonth) return TRUE;
    if (currmonth == compmonth) {
      if (currday > compday) return TRUE;
    }
  }

  return FALSE;
}

#include <ent2api.h>
extern REG CORE_GetREG(void);

static void SetNetIdent (void)

{
  static const char kRevision [] = "$Revision: " SEQ_APP_VER;
  unsigned int major = 0, minor = 0;
  char progname [80];
  char buf [128];
  CharPtr s, bf, pn;
  int n;
  int res;
  size_t len;

  pn = (CharPtr) progname;
  bf = (CharPtr) buf;
  StringCpy (pn, "Sequin");
  s = StringChr (kRevision, ':');
  if (s != 0) {
    res = sscanf (s + 1, "%u.%u%n", &major, &minor, &n);
    if (res >= 2 && n > 0) {
      len = StringLen (pn);
      sprintf (pn + len, "-%u.%u", major, minor);
      REG_Set (CORE_GetREG (), DEF_CONN_REG_SECTION,
               REG_CONN_ARGS, pn, eREG_Transient);
      pn [len] = '/';
    }
  }
  EntrezSetProgramName (pn);
  sprintf (bf, "User-Agent: %s\r\n", pn);
  REG_Set (CORE_GetREG (), DEF_CONN_REG_SECTION,
           REG_CONN_HTTP_USER_HEADER, bf, eREG_Transient);
}

extern CharPtr objPrtMemStr;

Int2 Main (void)

{
  AsnIoPtr       aip;
  BioseqPtr      bsp;
  BioseqSetPtr   bssp;
  Pointer        dataptr = NULL;
  BtnActnProc    fetchProc;
  CharPtr        filename = NULL;
  FILE           *fp;
  Int2           handled;
  Boolean        notaxid;
  SeqEntryPtr    oldsep;
  /*
  ObjMgrPtr      omp;
  OMProcControl  ompc;
  ObjMgrProcPtr  ompp;
  */
  OMUserDataPtr  omudp;
  PaneL          p;
  CharPtr        ptr;
  SeqEntryPtr    sep;
  Int4           smartPort = 0; 
  Char           str [80];
  Int2           val;
  WindoW         w;
#if defined(OS_UNIX) || defined(WIN_MOTIF)
  Int2           i;
#endif
#ifdef WIN_MOTIF
  RecT           r;
#endif
#if defined(OS_MAC) && !defined(OS_UNIX_DARWIN)
  long           sysVer;
#endif

  ErrSetFatalLevel (SEV_MAX);
  ErrClearOptFlags (EO_SHOW_USERSTR);
  ProcessUpdatesFirst (FALSE);

  SetNetIdent ();

  UseLocalAsnloadDataAndErrMsg ();
  ErrPathReset ();

#if defined(OS_MAC) && !defined(OS_UNIX_DARWIN)
  if ( Gestalt (gestaltSystemVersion, &sysVer) == noErr) {
    /* system version in low order word is hexadecimal */
    if (sysVer >= 4096) {
      Message (MSG_OK, "You are running on MacOS X and should use the native version of Sequin, not SequinOS9");
    }
  }
#endif

  if (MoreThanYearOld ()) {
    Message (MSG_OK, "This copy of Sequin is more than a year old.  Please download the current version.");
  }

  ReadSettings ();

  sprintf (str, "Sequin Application Version %s", SEQUIN_APPLICATION);
  SEQUIN_VERSION = StringSave (str);

  ptr = "Standard Release";
  if (useEntrez || useBlast) {
    ptr = "Network Aware";
  }
  if (indexerVersion) {
    ptr = "Indexer Services";
  }
  if (genomeCenter != NULL) {
    ptr = "Genome Center";
  }
  if (indexerVersion) {
    sprintf (str, "%s [%s - %s]", ptr, date_of_compilation, time_of_compilation);
  } else {
    sprintf (str, "%s [%s]", ptr, date_of_compilation);
  }
  SEQUIN_SERVICES = StringSave (str);

  if (! indexerVersion) {
    ErrSetLogfile (NULL, 0);
  }

  w = FixedWindow (-50, -33, -10, -10, "Sequin", NULL);
  p = SimplePanel (w, AboutBoxWidth (), AboutBoxHeight (), DrawAbout);
  Show (w);
#ifdef WIN_MOTIF
  Select (w);
  ObjectRect (p, &r);
  Select (p);
  InsetRect (&r, 3, 3);
  InvalRect (&r);
#endif
  Update ();

  workbenchMode = FALSE;
  subtoolMode = FALSE;
  stdinMode = FALSE;
  bioseqsetMode = FALSE;
  binseqentryMode = FALSE;
  entrezMode = FALSE;
  nohelpMode = FALSE;
  subtoolDatatype = 0;
  subtoolEntityID = 0;
  leaveAsOldAsn = FALSE;

#if defined(OS_UNIX) || defined(WIN_MOTIF)
  {{
    Nlm_Int4         argc = GetArgc();
    Nlm_CharPtr PNTR argv = GetArgv();
    for (i = 1;  i < argc;  i++)
      {
        if      (StringCmp (argv[i], "-s") == 0)
          subtoolMode = TRUE;
        else if (StringCmp (argv[i], "-w") == 0)
          workbenchMode = TRUE;
        else if (StringCmp (argv[i], "-x") == 0)
          stdinMode = TRUE;
        else if (StringCmp (argv[i], "-f") == 0 && i + 1 < argc) {
          stdinMode = TRUE;
          filename = argv [i + 1];
        } else if (StringCmp (argv[i], "-bse") == 0)
          binseqentryMode = TRUE;
        else if (StringCmp (argv[i], "-e") == 0)
          entrezMode = TRUE;
        else if (StringCmp (argv[i], "-h") == 0)
          nohelpMode = TRUE;
        /*
        else if (StringCmp (argv[i], "-udv") == 0)
          useUdv = TRUE;
        */
        else if (StringCmp (argv[i], "-oldgph") == 0)
          useOldGraphicView = TRUE;
        else if (StringCmp (argv[i], "-oldaln") == 0)
          useOldAlignmentView = TRUE;
        else if (StringCmp (argv[i], "-oldseq") == 0)
          useOldSequenceView = TRUE;
        else if (StringCmp (argv[i], "-gc") == 0) {
            indexerVersion = FALSE;
            extraServices = TRUE;
            genomeCenter = StringSave ("genome center tag goes here");
        }
        if (StringCmp (argv[i], "-b") == 0) {
          bioseqsetMode = TRUE;
        } else if (StringCmp (argv[i], "-oldaln") == 0) {
          newAlignReader = FALSE;
        } else if (StringCmp (argv[i], "-oldasn") == 0) {
          leaveAsOldAsn = TRUE;
        } else if (StringCmp (argv[i], "-oldsource") == 0) {
          SetAppProperty ("OldFlatfileSource", (void *) 1024);
        } else if (StringCmp (argv[i], "-a") == 0) {
          gphviewscorealigns = TRUE;
        } else if (StringCmp (argv[i], "-y") == 0) {
          backupMode = TRUE;
        } else if (StringCmp (argv[i], "-noseqfetch") == 0) {
          useSeqFetch = FALSE;
        } else if (StringCmp (argv[i], "-nolocalfetch") == 0) {
          useLocal = FALSE;
        } else if (StringCmp (argv[i], "-noseqidlookup") == 0) {
          useIdLookup = FALSE;
        }
#ifdef USE_SMARTNET
        else if (StringCmp (argv[i], "-ds") == 0) {
          dirsubMode = TRUE;
          nohelpMode = TRUE;
        }
        else if (StringNCmp (argv[i], "-z", 2) == 0) {
          smartnetMode = TRUE;
          dirsubMode = TRUE;
          if(*(argv[i]+2) != NULLB)
            smartPort = atoi(argv[i]+2);
          else
            smartPort = SM_SERVER_PORT;
        }
#endif
      }
  }}
#endif

#ifdef WIN_MSWIN
  {{
    Nlm_Int2         i;
    Nlm_Int4         argc = GetArgc();
    Nlm_CharPtr PNTR argv = GetArgv();
    for (i = 1;  i < argc;  i++)
      {
        if (StringCmp (argv[i], "-x") == 0)
          stdinMode = TRUE;
        else if (StringCmp (argv[i], "-f") == 0 && i + 1 < argc) {
          stdinMode = TRUE;
          filename = argv [i + 1];
        } else if (StringCmp (argv[i], "-e") == 0)
          entrezMode = TRUE;
        else if (StringCmp (argv[i], "-h") == 0)
          nohelpMode = TRUE;
        else if (StringCmp (argv[i], "-noseqfetch") == 0)
          useSeqFetch = FALSE;
        else if (StringCmp (argv[i], "-nolocalfetch") == 0)
          useLocal = FALSE;
        else if (StringCmp (argv[i], "-noseqidlookup") == 0)
          useIdLookup = FALSE;
#ifdef USE_SMARTNET
        else if (StringNCmp (argv[i], "-z", 2) == 0) {
          smartnetMode = TRUE;
          dirsubMode = TRUE;
          if(*(argv[i]+2) != NULLB)
            smartPort = atoi(argv[i]+2);
          else
            smartPort = SM_SERVER_PORT;
        }
#endif
      }
  }}
#endif

  if (subtoolMode || smartnetMode) {
    Char  tmp [PATH_MAX];
    ProgramPath (tmp, sizeof (tmp));
    ptr = StringRChr (tmp, DIRDELIMCHR);
    if (ptr != NULL) {
      ptr++;
    }
    if (smartPort > 0) {
      sprintf (str, "- %s %s (%ld) [%s - %s]", ptr, SEQUIN_APPLICATION, (long) smartPort, date_of_compilation, time_of_compilation);
    } else {
      sprintf (str, "- %s %s [%s - %s]", ptr, SEQUIN_APPLICATION, date_of_compilation, time_of_compilation);
    }
    SetAppProperty ("SmartSequinTimeStampTitle", StringSave (str));
  }

  if (GetSequinAppParam ("SETTINGS", "STDINMODE", NULL, str, sizeof (str))) {
    if (StringICmp (str, "TRUE") == 0) {
      stdinMode = TRUE;
    }
  }
  if (GetSequinAppParam ("SETTINGS", "HIDEHELPFORM", NULL, str, sizeof (str))) {
    if (StringICmp (str, "TRUE") == 0) {
      nohelpMode = TRUE;
    }
  }

  WatchCursor ();

  SetTitle (w, "Loading parse tables");
  if (! AllObjLoad ()) {
    ArrowCursor ();
    Message (MSG_FATAL, "AllObjLoad failed");
    return 0;
  }
  if (! SubmitAsnLoad ()) {
    ArrowCursor ();
    Message (MSG_FATAL, "SubmitAsnLoad failed");
    return 0;
  }
  if (! objprojAsnLoad ()) {
    ArrowCursor ();
    Message (MSG_FATAL, "ProjectAsnLoad failed");
    return 0;
  }

  SetTitle (w, "Loading print templates");
  /* objprt.prt still needed for Edit Citations button and Desktop view */
  if (! PrintTemplateSetLoadEx ("objprt.prt", objPrtMemStr)) {
    ArrowCursor ();
    Message (MSG_FATAL, "PrintTemplateSetLoad objprt.prt failed");
    return 0;
  }

  SetTitle (w, "Loading sequence alphabet converter");
  if (! SeqCodeSetLoad ()) {
    ArrowCursor ();
    Message (MSG_FATAL, "SeqCodeSetLoad failed");
    return 0;
  }

  SetTitle (w, "Loading organism table");
  if (! LoadOrganismTable ()) {
    ArrowCursor ();
    Message (MSG_POSTERR, "LoadOrganismTable failed");
  }

  SetTitle (w, "Loading genetic code table");
  if (! GeneticCodeTableLoad ()) {
    ArrowCursor ();
    Message (MSG_FATAL, "GeneticCodeTableLoad failed");
    return 0;
  }

  SetTitle (w, "Loading feature definitions");
  if (! FeatDefSetLoad ()) {
    ArrowCursor ();
    Message (MSG_FATAL, "FeatDefSetLoad failed");
    return 0;
  }

  SetupGeneticCodes ();

  GeneticCodeSingletonInit ();

  if (! SetupPrintOptions ()) {
    ArrowCursor ();
    Message (MSG_FATAL, "StdPrintOptionsNew failed");
    return 0;
  }

/*
#ifndef WIN16
  SetTitle (w, "Loading structure dictionary");
  if (OpenMMDBAPI ((POWER_VIEW ^ FETCH_ENTREZ), NULL)) {
    prgdDict = GetPRGDDictionary ();
    if (BiostrucAvail ()) {
      Cn3DWin_Entrez(NULL, useEntrez);
    }
  }
#endif
*/

  SetTitle (w, "Creating menus");
  SetupDesktop ();
  SetupCommonFonts ();

  InitSequinExtras ();

  VSMFileInit ();
  VSeqMgrInit (FALSE);
  WatchCursor ();

  /* register fetch functions */

  if (useEntrez) {
    /* EntrezBioseqFetchEnable ("Sequin", FALSE); */
    if (useSeqFetch) {
      /* ID1BioseqFetchEnable ("Sequin", FALSE); */
      PubSeqFetchEnable ();
      PubMedFetchEnable ();
    } else if (useIdLookup) {
      PubSeqFetchEnableEx (FALSE, TRUE, TRUE, TRUE, TRUE);
      PubMedFetchEnable ();
    }
  }

#ifdef USE_SMARTNET
  if (dirsubMode) {
    if (useSeqFetch) {
      /* DirSubFetchEnable (); */
      TPASmartFetchEnable ();
      SmartFetchEnable ();
      HUPFetchEnable ();
    }
  }
#endif

  if (useLocal) {
    LocalSeqFetchInit (FALSE);
  }

/*
#ifdef USE_SMARTNET
  if (dirsubMode) {
    if (only_use_smart) {
      SmartFetchEnable ();
    } else {
      DirSubFetchEnable ();
      SmartFetchEnable ();
      TPASmartFetchEnable ();
    }
  }
#endif

  if (! only_use_smart) {
    if (useEntrez) {
      if (useSeqFetch) {
        PubSeqFetchEnable ();
        PubMedFetchEnable ();
      } else {
        PubSeqFetchEnableEx (FALSE, TRUE, TRUE, TRUE, TRUE);
        PubMedFetchEnable ();
      }
    }

    if (useLocal) {
      LocalSeqFetchInit (FALSE);
    }
  }
*/

#ifdef WIN_MAC
  SetDeactivate (NULL, MacDeactProc);
  SetupMacMenus ();
#endif

  SetTitle (w, "Creating window");
  InitMuskStyles ();

  startupStyle = 0;
  if (GetSequinAppParam ("SETTINGS", "DEFAULTSTYLE", NULL, str, sizeof (str))) {
    if (StrToInt (str, &val) && val >= 0) {
      startupStyle = val;
      SetMuskCurrentSt (GetMuskStyleName (val));
    }
  }

  /* SequinCheckSocketsProc is called by SubtoolModeTimerProc,
  so that can safely override the Metronome call */

  Metronome (SequinCheckSocketsProc);

  subtoolEntityID = 0;
  subtoolTimerLimit = 100;
  subtoolTimerCount = 0;
  subtoolRecordDirty = FALSE;
  if (subtoolMode || smartnetMode) {
    if (GetSequinAppParam ("SETTINGS", "TIMERLIMIT", NULL, str, sizeof (str))) {
      if (StrToInt (str, &val) && val >= 0) {
        subtoolTimerLimit = val;
      }
    }
  }

/* -------------------------- SmartNet -------------------------- */
  if(smartnetMode) {
/*
    omp = ObjMgrGet();
    MemSet((Pointer)(&ompc), 0, sizeof(OMProcControl));
    ompc.input_entityID = 0;
    ompc.input_itemID = 0;
    ompc.input_itemtype = OBJ_BIOSEQ;
    ompc.input_data = NULL;
    ompp = NULL;
    while (procval != OM_MSG_RET_DONE &&
           (ompp = ObjMgrProcFindNext (omp, OMPROC_VIEW, OBJ_BIOSEQ, OBJ_BIOSEQ, ompp)) != NULL) {
      if (! ompp->subinputtype) {
        ompc.proc = ompp;
        procval = (*(ompp->func))((Pointer)&ompc);
        if (procval == OM_MSG_RET_ERROR) {
          ErrShow ();
        }
      }
    }
*/

#ifdef USE_SMARTNET
    Remove (w);
    ArrowCursor ();

    if(!SMListenRequests(NULL,
                         SMReadBioseqObj,
                         smartPort)) {
      Message(MSG_OK, "Cannot start accept thread");
      QuitProgram();
    }
    subtoolTimerCount = 0;
    subtoolRecordDirty = FALSE;
    Metronome (SubtoolModeTimerProc);
    ProcessEvents();

    /* No more code in Main() for smartnet Mode */
    return 0;

#endif
/* --------------------------------------------------------------- */

  } else if (subtoolMode || stdinMode || binseqentryMode) {
    Remove (w);
    ArrowCursor ();
    if (subtoolMode) {
      dataptr = SubtoolModeAsnTextFileRead ("stdin", &subtoolDatatype, &subtoolEntityID);
    } else if (binseqentryMode) {
      dataptr = NULL;
      aip = AsnIoOpen ("stdin", "rb");
      SeqMgrHoldIndexing (TRUE);
      sep = SeqEntryAsnRead (aip, NULL);
      SeqMgrHoldIndexing (FALSE);
      AsnIoClose (aip);
      if (sep != NULL) {
        if (sep->choice == 1) {
          subtoolDatatype = OBJ_BIOSEQ;
          dataptr = (Pointer) sep->data.ptrvalue;
        } else if (sep->choice == 2) {
          subtoolDatatype = OBJ_BIOSEQSET;
          dataptr = (Pointer) sep->data.ptrvalue;
        }
        if (dataptr != NULL) {
          subtoolEntityID = ObjMgrRegister (subtoolDatatype, dataptr);
        }
      }
    } else {
      if (! StringHasNoText (filename)) {
        fp = FileOpen (filename, "r");
      } else {
        fp = FileOpen ("stdin", "r");
      }
      dataptr = ReadAsnFastaOrFlatFile (fp, &subtoolDatatype,  &subtoolEntityID, FALSE, FALSE, TRUE, FALSE);
      FileClose (fp);
      /*
      dataptr = ObjMgrGenericAsnTextFileRead ("stdin", &subtoolDatatype, &subtoolEntityID);
      */
      if (dataptr == NULL) {
        fseek (stdin, 0, SEEK_SET);
        aip = AsnIoOpen ("stdin", "rb");
        SeqMgrHoldIndexing (TRUE);
        sep = SeqEntryAsnRead (aip, NULL);
        SeqMgrHoldIndexing (FALSE);
        AsnIoClose (aip);
        if (sep != NULL) {
          if (sep->choice == 1) {
            subtoolDatatype = OBJ_BIOSEQ;
            dataptr = (Pointer) sep->data.ptrvalue;
          } else if (sep->choice == 2) {
            subtoolDatatype = OBJ_BIOSEQSET;
            dataptr = (Pointer) sep->data.ptrvalue;
          }
          if (dataptr != NULL) {
            subtoolEntityID = ObjMgrRegister (subtoolDatatype, dataptr);
          }
        }
      }
    }
    if (dataptr != NULL && subtoolEntityID > 0) {
      if (subtoolDatatype == OBJ_SEQSUB || subtoolDatatype == OBJ_SEQENTRY ||
          subtoolDatatype == OBJ_BIOSEQ || subtoolDatatype == OBJ_BIOSEQSET) {
        WatchCursor ();
        sep = GetTopSeqEntryForEntityID (subtoolEntityID);
        if (sep == NULL) {
          sep = SeqEntryNew ();
          if (sep != NULL) {
            if (subtoolDatatype == OBJ_BIOSEQ) {
              bsp = (BioseqPtr) dataptr;
              sep->choice = 1;
              sep->data.ptrvalue = bsp;
              SeqMgrSeqEntry (SM_BIOSEQ, (Pointer) bsp, sep);
            } else if (subtoolDatatype == OBJ_BIOSEQSET) {
              bssp = (BioseqSetPtr) dataptr;
              sep->choice = 2;
              sep->data.ptrvalue = bssp;
              SeqMgrSeqEntry (SM_BIOSEQSET, (Pointer) bssp, sep);
            } else {
              sep = SeqEntryFree (sep);
            }
          }
          sep = GetTopSeqEntryForEntityID (subtoolEntityID);
        }
        if (sep != NULL) {
          if (EntrezASN1Detected (sep)) {
            ErrPostEx (SEV_WARNING, 0, 0, "This record was retrieved from Entrez");
          }
          if (! leaveAsOldAsn) {
            if (subtoolMode) {
              MySeqEntryToAsn3 (sep, TRUE, FALSE, TRUE);

              /* now instantiating protein titles */
              InstantiateProteinTitles (subtoolEntityID, NULL);
            } else {
              notaxid = FALSE;
              VisitBioSourcesInSep (sep, (Pointer) &notaxid, LookForTaxonID);
              if (notaxid) {
                MySeqEntryToAsn3 (sep, TRUE, FALSE, FALSE);
              }
            }
          }
          if (subtoolMode) {
            if (sep->choice == 2 && sep->data.ptrvalue != NULL) {
              bssp = (BioseqSetPtr) sep->data.ptrvalue;
              if (bssp->_class != 7) {
                Message (MSG_OK, "WARNING: Converting from bioseq set of class %d", (int) bssp->_class);
              }
              bssp->_class = 7;
            }
            if (FileLength (SEQUIN_EDIT_BACK_FILE) > 0) {
              if (Message (MSG_YN, "Restore from backup?") == ANS_YES) {
                oldsep = RestoreFromFile (SEQUIN_EDIT_BACK_FILE);
                ReplaceSeqEntryWithSeqEntry (sep, oldsep, TRUE);
                subtoolEntityID = ObjMgrGetEntityIDForChoice (sep);
              }
            }
          }
        }
        if (stdinMode) {
          seqviewprocs.filepath = filename;
        }
        seqviewprocs.forceSeparateViewer = TRUE;
        SeqEntrySetScope (NULL);
        handled = GatherProcLaunch (OMPROC_VIEW, FALSE, subtoolEntityID, 1,
                                    OBJ_BIOSEQ, 0, 0, OBJ_BIOSEQ, 0);
        seqviewprocs.filepath = NULL;
        ArrowCursor ();
        subtoolTimerCount = 0;
        subtoolRecordDirty = FALSE;
        if (handled != OM_MSG_RET_DONE || handled == OM_MSG_RET_NOPROC) {
          Message (MSG_FATAL, "Unable to launch viewer.");
          CleanupSequin ();
          return 0;
        } else {
          if (! nohelpMode) {
            SetTitle (w, "Creating help window");
            if (helpForm == NULL) {
              helpForm = CreateHelpForm (-95, -5, "Sequin Help", "sequin.hlp",
                                         HideHelpForm, HelpActivateProc);
            }
          }
          SendHelpScrollMessage (helpForm, "Editing the Record", NULL);
        }
        ObjMgrSetOptions (OM_OPT_FREE_IF_NO_VIEW, subtoolEntityID);
        if (! subtoolMode) {
          ObjMgrSetDirtyFlag (subtoolEntityID, TRUE);
        }
        if (subtoolMode) {
          Metronome (SubtoolModeTimerProc);
          omudp = ObjMgrAddUserData (subtoolEntityID, 0, 0, 0);
          if (omudp != NULL) {
            omudp->messagefunc = SubtoolModeMsgFunc;
          }
          subtoolRecordDirty = TRUE;
        }
      } else {
        Message (MSG_FATAL, "Unable to recognize ASN.1 type %d", (int) subtoolDatatype);
        CleanupSequin ();
        return 0;
      }
    } else {
      Message (MSG_FATAL, "Unable to read ASN.1 file");
      CleanupSequin ();
      return 0;
    }
  } else if (workbenchMode) {
    Remove (w);
    ArrowCursor ();
  } else if (entrezMode) {
    Remove (w);
    ArrowCursor ();
  } else {
    if (! nohelpMode) {
      SetTitle (w, "Creating help window");
      if (helpForm == NULL) {
        helpForm = CreateHelpForm (-95, -5, "Sequin Help", "sequin.hlp",
                                   HideHelpForm, HelpActivateProc);
      }
    }
    SetTitle (w, "Creating initial window");
    fetchProc = NULL;
    if (allowDownload) {
      fetchProc = FetchFromNet;
    }
    if (genomeCenter != NULL) {
      startupForm = CreateStartupForm (-5, -67, "Welcome to Sequin",
                                       StartFa2htgs, StartPhrap, BuildContig,
                                       StartNew, ReadOld, fetchProc,
                                       ShowHelp, CreateSubmissionTemplate, DoQuit, StartupActivateProc);
    } else {
      startupForm = CreateStartupForm (-5, -67, "Welcome to Sequin",
                                       NULL, NULL, NULL, StartNew, ReadOld, fetchProc,
                                       ShowHelp, NULL, DoQuit, StartupActivateProc);
    }
    globalFormatBlock.seqPackage = SEQ_PKG_SINGLE;
    globalFormatBlock.seqFormat = SEQ_FMT_FASTA;
    globalFormatBlock.numSeqs = 0;
    globalFormatBlock.submType = SEQ_ORIG_SUBMISSION;
    Remove (w);
    ArrowCursor ();
  }

  if (backupMode) {
    Metronome (BackupModeTimerProc);
  }

  if (subtoolMode || stdinMode || binseqentryMode) {
  } else if (workbenchMode) {
  } else if (entrezMode) {
    /*
    MakeTermListForm ();
    if (termListForm != NULL) {
      Show (termListForm);
      Select (termListForm);
      Update ();
      MakeDocSumForm ();
      if (docSumForm != NULL) {
      } else {
        Message (MSG_FATAL, "Unable to create document window");
        CleanupSequin ();
        return 0;
      }
    } else {
      Message (MSG_FATAL, "Unable to create term list window");
      CleanupSequin ();
      return 0;
    }
    */
    Message (MSG_FATAL, "This mode is obsolete");
    CleanupSequin ();
    return 0;
  } else if (startupForm != NULL) {
    Show (startupForm);
    Select (startupForm);
    Update ();
    initSubmitForm = CreateInitSubmitterForm (-5, -67, "Submitting Authors",
                                              GetFormat, BackToStartup,
                                              SubmitBlockActivateProc);
    formatForm = CreateFormatForm (-5, -67, "Sequence Format",
                                   GetOrgAndSeq, BackToSubmitter, FormatActivateProc);
    if (! nohelpMode) {
      Update ();
      Show (helpForm);
      Select (helpForm);
      Update ();
      SendHelpScrollMessage (helpForm, "Introduction", NULL);
    }
  } else {
    Message (MSG_FATAL, "Unable to create window.");
    CleanupSequin ();
    return 0;
  }

#if defined(WIN_MAC)  ||  defined(WIN_MSWIN)
  RegisterDropProc (SequinOpenMimeFile);
  RegisterResultProc (SequinOpenResultFile);
#endif

  if (workbenchMode) {
    VSeqMgrRun ("NCBI Desktop", "Sequin Workbench");
  } else {
    ProcessEvents ();
  }

  WatchCursor ();
  val = GetMuskCurrentSt ();
  if (val < 0) {
    val = 0;
  }
  if (val >= GetMuskTotalSt ()) {
    val = 0;
  }
  if (val != startupStyle) {
    if (val > 0) {
      sprintf (str, "%d", (int) val);
      SetSequinAppParam ("SETTINGS", "DEFAULTSTYLE", str);
    } else {
      SetSequinAppParam ("SETTINGS", "DEFAULTSTYLE", "0");
    }
  }

  Remove (startupForm);
  Remove (initSubmitForm);
  Remove (formatForm);
  Remove (helpForm);

  FreeSqnTempFiles ();

  ArrowCursor ();
  Update ();

  CleanupSequin ();

  return 0;
}

extern SubmitBlockPtr ConvertSequinBlockToSubmitBlock (SequinBlockPtr sqp);

extern Boolean ExportSubmitterBlockTemplate (SeqEntryPtr sep, SeqDescrPtr sdp)
{
  SeqSubmitPtr ssp;
  Char         path [PATH_MAX];
  AsnIoPtr     aip;
  FILE         *fp;
  
  if (globalsbp == NULL)
  {
    Message (MSG_ERROR, "No submit block!");
    return FALSE;
  }
  
  if (!GetOutputFileName (path, sizeof (path), NULL))
  {
    return FALSE;
  }
  
  fp = FileOpen (path, "w");
  if (fp == NULL)
  {
    Message (MSG_ERROR, "Unable to open %s", path);
    return FALSE;
  }
  
  ssp = SeqSubmitNew ();
  if (ssp == NULL)
  {
    Message (MSG_ERROR, "Unable to allocate memory for seq-submit");
    return FALSE;
  }

  ssp->datatype = 1;
  ssp->data = (Pointer) sep;

  /* need to store copy of globalbsp in form */
  PointerToForm (initSubmitForm, globalsbp);
  ssp->sub = ConvertSequinBlockToSubmitBlock (globalsbp);
  /* get globalbsp back */
  globalsbp = FormToPointer (initSubmitForm);
  aip = AsnIoNew(ASNIO_TEXT_OUT, fp, NULL, NULL, NULL);
  SeqSubmitAsnWrite(ssp, aip, NULL);
  
	AsnIoFlush(aip);
  AsnIoReset(aip);
  if (sdp != NULL)
  {
    SeqDescAsnWrite (sdp, aip, NULL);  
	  AsnIoFlush(aip);
    AsnIoReset(aip);
  }
  AsnIoClose (aip);
  
  ssp = SeqSubmitFree (ssp);

  globalsbp = SequinBlockFree (globalsbp);
  Hide (initSubmitForm);
  Hide (formatForm);
  Update ();
  Show (startupForm);
  Select (startupForm);
  SendHelpScrollMessage (helpForm, "Introduction", NULL);
  Update ();    
  return TRUE;
}


extern void MakeBadSpecificHostValueTable (IteM i)
{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;
  LogInfoPtr   lip;
  ValNodePtr    misspelled = NULL, bad_caps = NULL, ambiguous = NULL, unrecognized = NULL;

#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;


  Taxon3ValidateSpecificHostsInSeqEntry (sep, &misspelled, &bad_caps, &ambiguous, &unrecognized);

  if (misspelled == NULL && bad_caps == NULL && ambiguous == NULL && unrecognized == NULL)
  {
    Message (MSG_OK, "No bad specific-host values found!");
    return;
  }
  lip = OpenLog ("Bad Specific-Host Values");

  if (misspelled != NULL) {
    fprintf (lip->fp, "Mis-spelled Specific-Host Values\n");
    lip->data_in_log = WriteBadSpecificHostTable (misspelled, lip->fp);
    misspelled = ValNodeFree (misspelled);
  }
  if (bad_caps != NULL) {
    fprintf (lip->fp, "Incorrectly Capitalized Specific-Host Values\n");
    lip->data_in_log = WriteBadSpecificHostTable (bad_caps, lip->fp);
    bad_caps = ValNodeFree (bad_caps);
  }
  if (ambiguous != NULL) {
    fprintf (lip->fp, "Ambiguous Specific-Host Values\n");
    lip->data_in_log = WriteBadSpecificHostTable (ambiguous, lip->fp);
    ambiguous = ValNodeFree (ambiguous);
  }
  if (unrecognized != NULL) {
    fprintf (lip->fp, "Unrecognized Specific-Host Values\n");
    lip->data_in_log = WriteBadSpecificHostTable (unrecognized, lip->fp);
    unrecognized = ValNodeFree (unrecognized);
  }
  CloseLog (lip);
  lip = FreeLog (lip);
}


typedef struct updatefeaturesform {
  FORM_MESSAGE_BLOCK
  DialoG new_features;
  DialoG old_features;
} UpdateFeaturesFormData, PNTR UpdateFeaturesFormPtr;


static void MoveFeatureToReplacement (ButtoN b)
{
  UpdateFeaturesFormPtr f;
  ValNodePtr            features;

  f = (UpdateFeaturesFormPtr) GetObjectExtra (b);
  if (f == NULL) return;

  features = RemoveSelectedFeaturesFromList (f->old_features);
  if (!AddFeaturesToReplaceList(f->new_features, features)) {
    features = ValNodeFree (features);
  }
}


static void MoveAllFeaturesToReplacement (ButtoN b)
{
  UpdateFeaturesFormPtr f;
  ValNodePtr            features;

  f = (UpdateFeaturesFormPtr) GetObjectExtra (b);
  if (f == NULL) return;

  features = (ValNodePtr) DialogToPointer (f->old_features);
  PointerToDialog (f->old_features, NULL);
  if (!AddFeaturesToReplaceList(f->new_features, features)) {
    features = ValNodeFree (features);
  }
}


static void AutomatchFeaturesForReplacement (ButtoN b) 
{
  UpdateFeaturesFormPtr f;
  ValNodePtr            existing_features;

  f = (UpdateFeaturesFormPtr) GetObjectExtra (b);
  if (f == NULL) return;

  existing_features = DialogToPointer (f->old_features);
  PointerToDialog (f->old_features, NULL);
  if (AutomatchFeatures (f->new_features, &existing_features)) {
    PointerToDialog (f->old_features, existing_features);
  }

}


static void MoveFeatureToSelection (ButtoN b)
{
  UpdateFeaturesFormPtr f;
  ValNodePtr            features;

  f = (UpdateFeaturesFormPtr) GetObjectExtra (b);
  if (f == NULL) return;

  features = RemoveFeaturesFromReplaceList (f->new_features);
  AddFeaturesToList (f->old_features, features);
}


static void DoUpdateFeatures (ButtoN b)
{
  UpdateFeaturesFormPtr f;
  ValNodePtr            feature_list;

  f = (UpdateFeaturesFormPtr) GetObjectExtra (b);
  if (f == NULL) return;

  feature_list = DialogToPointer (f->new_features);
  ActOnFeatureReplaceList (feature_list);
  feature_list = FeatureReplaceListFree (feature_list);
  DeleteMarkedObjects (f->input_entityID, 0, NULL);
  ObjMgrSetDirtyFlag (f->input_entityID, TRUE);
  ObjMgrSendMsg (OM_MSG_UPDATE, f->input_entityID, 0, 0);

  Remove (f->form);
}



static void ScrollToReplacementFeature (Pointer data)
{
  UpdateFeaturesFormPtr frm;
  SeqFeatPtr new_feat;

  frm = (UpdateFeaturesFormPtr) data;
  if (frm == NULL) {
    return;
  }

  /* find currently highlighted feature */
  new_feat = GetSelectedNewFeature (frm->new_features);

  /* find "auto" replacement feature */
  /* if found, highlight and scroll to it */
  /* otherwise scroll to first feature that overlaps this location */
  ScrollToMatchingFeatures (frm->old_features, new_feat);
}


extern void UpdateFeatures (IteM i)
{
  BaseFormPtr  bfp;
  SeqEntryPtr  sep;
  Pointer        dataptr;
  Uint2          datatype;
  FILE           *fp;
  Char           path [PATH_MAX];
  SeqAnnotPtr    sap;
  SeqFeatPtr     sfp;
  ValNodePtr     new_feat_list = NULL, no_bsp_list = NULL, old_item_list = NULL;
  Int4           leftmost = -1, rightmost = -1, new_left, new_right, tmp;
  BioseqPtr      bsp, last_bsp = NULL;
  SeqLocPtr        slp;
  WindoW           w;
  GrouP            h, g2, k, c;
  ButtoN           b;
  UpdateFeaturesFormPtr f;


#ifdef WIN_MAC
  bfp = currentFormDataPtr;
#else
  bfp = GetObjectExtra (i);
#endif
  if (bfp == NULL) return;
  sep = GetTopSeqEntryForEntityID (bfp->input_entityID);
  if (sep == NULL) return;

  path [0] = '\0';
  if (!GetInputFileName (path, sizeof (path), "", "TEXT")) {
    return;
  }

  fp = FileOpen (path, "r");
  if (fp == NULL) {
    Message (MSG_ERROR, "Unable to open file!");
    return;
  }
  dataptr = ReadAsnFastaOrFlatFile (fp, &datatype, NULL, FALSE, FALSE,
                                    TRUE, TRUE);
  FileClose (fp);
  if (dataptr == NULL || datatype != OBJ_SEQANNOT) {
    Message (MSG_ERROR, "File does not contain feature table!");
    return;
  }

  sap = (SeqAnnotPtr) dataptr;
  if (sap->type != 1) {
    Message (MSG_ERROR, "File does not contain feature table!");
    sap = SeqAnnotFree (sap);
    return;
  }

  SortSeqFeatInAnnot (sap);

  /* get intervals for features for each bioseq */
  for (sfp = sap->data; sfp != NULL; sfp = sfp->next) {
    bsp = BioseqFindFromSeqLoc (sfp->location);
    if (bsp == NULL) {    
      ValNodeAddPointer (&no_bsp_list, OBJ_SEQFEAT, sfp);
    } else {
      if (bsp != last_bsp) {   
        if (last_bsp != NULL) {
          slp = SeqLocIntNew (leftmost, rightmost, Seq_strand_plus, SeqIdDup (SeqIdFindWorst (bsp->id)));
          ValNodeLink (&old_item_list, ListFeaturesOverlappingLocation (last_bsp, slp, 0, 0));
          slp = SeqLocFree (slp);
        }   
        last_bsp = bsp;
        leftmost = -1;
        rightmost = -1;
      }
      new_left = SeqLocStart (sfp->location);
      new_right = SeqLocStop (sfp->location);
      if (new_left > new_right) {
        tmp = new_left;
        new_left = new_right;
        new_right = tmp;
      }
      if (leftmost == -1 || new_left < leftmost) {
        leftmost = new_left;
      }
      if (rightmost == -1 || new_right > rightmost) {
        rightmost = new_right;
      }
    }
  }
  if (last_bsp != NULL) {
    slp = SeqLocIntNew (leftmost, rightmost, Seq_strand_plus, SeqIdDup (SeqIdFindWorst (bsp->id)));
    ValNodeLink (&old_item_list, ListFeaturesOverlappingLocation (last_bsp, slp, 0, 0));
    slp = SeqLocFree (slp);
  }

  if (no_bsp_list != NULL) {
    Message (MSG_ERROR, "%d features in table are not found on a Bioseq in this record!", ValNodeLen (no_bsp_list));
    no_bsp_list = ValNodeFree (no_bsp_list);
  }

  new_feat_list = FeatureReplaceListFromSeqAnnot (sap);
  sap = SeqAnnotFree (sap);
 
  if (new_feat_list == NULL) {
    Message (MSG_ERROR, "No features found!");
    return;
  }   

  /* Now create dialog to allow user to select new features to import and existing features to delete */
  f = (UpdateFeaturesFormPtr) MemNew (sizeof (UpdateFeaturesFormData));
  w = FixedWindow (-50, -33, -10, -10, "Update Features", StdCloseWindowProc);
  SetObjectExtra (w, f, StdCleanupFormProc);
  f->form = (ForM) w;
  f->input_entityID = bfp->input_entityID;
  h = HiddenGroup (w, -1, 0, NULL);
  SetGroupSpacing (h, 10, 10);

  g2 = HiddenGroup (h, 3, 0, NULL);
  StaticPrompt (g2, "Features to Import", 0, 0, programFont, 'c');
  StaticPrompt (g2, "", 0, 0, programFont, 'l');
  StaticPrompt (g2, "Existing Features that Overlap this Interval", 0, 0, programFont, 'c');
  f->new_features = FeatureReplaceListDialog (g2, 400, ScrollToReplacementFeature, f);
  PointerToDialog (f->new_features, new_feat_list);
  new_feat_list = FeatureReplaceListFree (new_feat_list);
  k = HiddenGroup (g2, 0, 3, NULL);
  b = PushButton (k, "<=", MoveFeatureToReplacement);
  SetObjectExtra (b, f, NULL);
  b = PushButton (k, "<<=", MoveAllFeaturesToReplacement);
  SetObjectExtra (b, f, NULL);
  b = PushButton (k, "=>", MoveFeatureToSelection);
  SetObjectExtra (b, f, NULL);
  f->old_features = FeatureSelectListDialog (g2, 400);
  PointerToDialog (f->old_features, old_item_list);

  c = HiddenGroup (h, 4, 0, NULL);
  b = PushButton (c, "Accept", DoUpdateFeatures);
  SetObjectExtra (b, f, NULL);
  b = PushButton (c, "Automatch", AutomatchFeaturesForReplacement);
  SetObjectExtra (b, f, NULL);
  PushButton (c, "Cancel", StdCancelButtonProc);
  
  AlignObjects (ALIGN_CENTER, (HANDLE) g2, (HANDLE) c, NULL);
  Show (w);
  Select (w);
}
