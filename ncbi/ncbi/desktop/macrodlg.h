/*   macrodlg.h
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
* File Name:  macrodlg.h
*
* Author:  Colleen Bollin
*
* Version Creation Date:   11/23/2007
*
* $Revision: 1.20 $
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

#ifndef MACRODLG_H
#define MACRODLG_H

#include <objmacro.h>

NLM_EXTERN void LaunchMacroEditor (IteM i);
NLM_EXTERN void LaunchMacroEditorBaseForm (BaseFormPtr bfp);

NLM_EXTERN DialoG TabColumnConfigDialog 
(GrouP                    h,
 CharPtr                  title,
 Int4                     num_blank,
 Nlm_ChangeNotifyProc     change_notify,
 Pointer                  change_userdata);
NLM_EXTERN DialoG TabColumnConfigListDialog (GrouP h, ValNodePtr first_values, ValNodePtr blank_list, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata);
NLM_EXTERN void ChangeDataForTabColumnConfigListDialog (DialoG d, ValNodePtr first_values, ValNodePtr blank_list);
NLM_EXTERN DialoG MatchTypeDialog (GrouP g, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata);

NLM_EXTERN DialoG FeatureTypeDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata);
NLM_EXTERN DialoG StringConstraintDialog (GrouP h, CharPtr label, Boolean clear_btn, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata);
NLM_EXTERN DialoG ComplexConstraintDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata);
NLM_EXTERN void ChangeComplexConstraintFieldType (DialoG d, Uint2 qual_type, ValNodePtr rna_type, Int2 feat_type);
NLM_EXTERN DialoG ConstraintSetDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata);
NLM_EXTERN DialoG MolInfoBlockDialog (GrouP h, Boolean edit, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata);
NLM_EXTERN void SingleAECRMacroAction (Uint2 entityID, Boolean indexer_version, Uint1 AECR_action_type, Uint1 AECR_qual_type);
NLM_EXTERN void MacroApplyKeyword (Uint2 entityID, Boolean indexer_version);

NLM_EXTERN Uint2 TwoStepExistingText (Int4 num_found, Boolean non_text, Boolean allow_multi);

NLM_EXTERN ForM SingleParseAction (Uint2 entityID);

NLM_EXTERN DialoG MolinfoFieldChoiceDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata);

/* for adding features */
typedef struct applyfeaturedetails {
  Boolean add_mrna;
  ValNodePtr fields;
  ValNodePtr src_fields;
} ApplyFeatureDetailsData, PNTR ApplyFeatureDetailsPtr;

NLM_EXTERN DialoG ApplyFeatureDetailsDialog (GrouP h, Uint1 featdef_type, ApplyFeatureDetailsPtr details, Boolean indexer_version, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata);
NLM_EXTERN ApplyFeatureDetailsPtr ApplyFeatureDetailsNew (ApplyFeatureActionPtr action);
NLM_EXTERN ApplyFeatureDetailsPtr ApplyFeatureDetailsFree (ApplyFeatureDetailsPtr details); 
NLM_EXTERN DialoG ParseDstDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata);
NLM_EXTERN DialoG CapChangeDialog (GrouP h, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata);
NLM_EXTERN DialoG TextPortionDialog (GrouP h, Boolean inside, Nlm_ChangeNotifyProc change_notify, Pointer change_userdata);


#endif
