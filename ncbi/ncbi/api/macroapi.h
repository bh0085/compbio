/*   macro_i.h
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
* File Name:  macro_i.h
*
* Author:  Colleen Bollin
*
* Version Creation Date:   11/15/2007
*
* $Revision: 1.90 $
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

#ifndef _macroapi_h_
#define _macroapi_h_

#ifdef __cplusplus
extern "C" {
#endif

/* 
 * Some batch operations will be faster if information about the entire record is collected once
 * and reused.  The BatchExtra structure is where such data belongs.
 */
typedef struct batchextra {
  ValNodePtr cit_list; /* this contains a list of minimized pubs and the numbers to be used for citations referring
                        * to these pubs on a given Bioseq.
                        * If needed, it should be generated by GetCitListsForSeqEntry.
                        */
} BatchExtraData, PNTR BatchExtraPtr;

NLM_EXTERN BatchExtraPtr BatchExtraNew ();
NLM_EXTERN BatchExtraPtr BatchExtraFree (BatchExtraPtr b);

NLM_EXTERN ValNodePtr GetCitListsForSeqEntry (SeqEntryPtr sep);
NLM_EXTERN ValNodePtr PubSerialNumberListFree (ValNodePtr vnp);
/* GetCitationNumberForMinPub can be used to calculate the citation number to be used
 * for a minimized pub from a SeqFeat->cit on a given Bioseq.  pub_list should have
 * been generated by GetCitListsForSeqEntry and should be freed by PubSerialNumberListFree.
 */
NLM_EXTERN Int4 GetCitationNumberForMinPub (BioseqPtr bsp, ValNodePtr min_pub, ValNodePtr pub_list);
NLM_EXTERN ValNodePtr GetMinPubForCitationNumber (BioseqPtr bsp, Int4 number, ValNodePtr pub_list);

NLM_EXTERN FeatureFieldPtr FeatureFieldCopy (FeatureFieldPtr orig);
NLM_EXTERN FieldTypePtr FieldTypeCopy (FieldTypePtr orig);

NLM_EXTERN FieldTypePtr FieldTypeFromString (CharPtr str);
NLM_EXTERN Int4 GetFeatdefFromFeatureType (Int4 feature_type);
NLM_EXTERN Int4 GetFeatureTypeFromFeatdef (Int4 featdef);
NLM_EXTERN CharPtr GetFeatureNameFromFeatureType (Int4 feature_type);
NLM_EXTERN Int4 GetFeatureTypeByName (CharPtr feat_name);
NLM_EXTERN void AddImportFeaturesToChoiceList (ValNodePtr PNTR feature_type_list);
NLM_EXTERN void AddAllFeaturesToChoiceList (ValNodePtr PNTR feature_type_list);
NLM_EXTERN CharPtr GetFeatQualName (Int4 featqual); 
NLM_EXTERN Int4 GetFeatQualByName (CharPtr qualname); 
NLM_EXTERN Int4 GetNumFeatQual (void);
NLM_EXTERN void AddAllFeatureFieldsToChoiceList (ValNodePtr PNTR field_list);
NLM_EXTERN CharPtr SummarizeFeatQual (ValNodePtr qual);
NLM_EXTERN CharPtr GetSourceQualName (Int4 srcqual);
NLM_EXTERN Int4 GetSourceQualTypeByName (CharPtr qualname);
NLM_EXTERN Int4 GetSrcQualFromSubSrcOrOrgMod (Int4 qual, Boolean is_org_mod);
NLM_EXTERN Int4 GetOrgModQualFromSrcQual (Int4 srcqual, Int4Ptr subfield);
NLM_EXTERN ValNodePtr GetSourceQualList (Boolean for_remove);
NLM_EXTERN Boolean IsNonTextSourceQual (Int4 srcqual);
NLM_EXTERN Boolean IsNonTextFieldType (FieldTypePtr field);
NLM_EXTERN TextFsaPtr GetOrgModSearch (void);
NLM_EXTERN Int4 GenomeFromSrcLoc (Int4 srcloc);
NLM_EXTERN Int4 SrcLocFromGenome (Int4 genome);
NLM_EXTERN CharPtr LocNameFromGenome (Int4 genome);
NLM_EXTERN Int4 GenomeFromLocName (CharPtr loc_name);
NLM_EXTERN ValNodePtr GetLocationList (Boolean for_remove); 
NLM_EXTERN Int4 OriginFromSrcOrig (Int4 srcorig);
NLM_EXTERN Int4 SrcOrigFromOrigin (Int4 origin);
NLM_EXTERN CharPtr OriginNameFromOrigin (Int4 origin);
NLM_EXTERN ValNodePtr GetOriginList (Boolean for_remove);
NLM_EXTERN BioSourcePtr GetBioSourceFromObject (Uint1 choice, Pointer data);
NLM_EXTERN ValNodePtr GetSourceQualSampleFieldList (SeqEntryPtr sep);
NLM_EXTERN ValNodePtr GetSourceQualSampleFieldListForSeqEntryList (ValNodePtr list);
NLM_EXTERN CharPtr CDSGeneProtNameFromField (Int4 field); 
NLM_EXTERN CharPtr CDSGeneProtFeatureNameFromFeatureType (Int4 feature_type);
NLM_EXTERN void AddAllCDSGeneProtFieldsToChoiceList (ValNodePtr PNTR field_list);
NLM_EXTERN void AddAllCDSGeneProtFeaturesToChoiceList (ValNodePtr PNTR field_list);
NLM_EXTERN FeatureFieldPtr FeatureFieldFromCDSGeneProtField (Uint2 cds_gene_prot_field);

NLM_EXTERN CharPtr BiomolNameFromBiomol (Int4 biomol);
NLM_EXTERN Int4 BiomolFromMoleculeType (Int4 molecule_type);
NLM_EXTERN ValNodePtr GetMoleculeTypeList (void);
NLM_EXTERN CharPtr TechNameFromTech (Int4 tech);
NLM_EXTERN Int4 TechFromTechniqueType (Int4 technique_type);
NLM_EXTERN ValNodePtr GetTechniqueTypeList (void);
NLM_EXTERN Int4 CompletenessFromCompletednessType (Int4 completedness_type);
NLM_EXTERN CharPtr CompletenessNameFromCompleteness (Int4 completeness); 
NLM_EXTERN ValNodePtr GetCompletednessTypeList (void);
NLM_EXTERN Int4 MolFromMoleculeClassType (Int4 moleculeclass_type);
NLM_EXTERN CharPtr MolNameFromMol (Int4 mol); 
NLM_EXTERN ValNodePtr GetMoleculeClassTypeList (void);
NLM_EXTERN Int4 TopologyFromTopologyType (Int4 topology_type);
NLM_EXTERN CharPtr TopologyNameFromTopology (Int4 topology);
NLM_EXTERN ValNodePtr GetTopologyTypeList (void);
NLM_EXTERN Int4 StrandFromStrandType (Int4 strand_type);
NLM_EXTERN CharPtr StrandNameFromStrand (Int4 strand);
NLM_EXTERN ValNodePtr GetStrandTypeList (void);
NLM_EXTERN Int4 Asn1BondTypeFromMacroBondType (Int4 macro_bond_type);
NLM_EXTERN Int4 MacroBondTypeFromAsn1BondType (Int4 asn1_bond_type); 
NLM_EXTERN CharPtr GetMacroBondTypeName (Int4 macro_bond_type);
NLM_EXTERN ValNodePtr GetBondTypeList (void);
NLM_EXTERN Int4 Asn1SiteTypeFromMacroSiteType (Int4 macro_site_type);
NLM_EXTERN Int4 MacroSiteTypeFromAsn1SiteType (Int4 asn1_site_type); 
NLM_EXTERN CharPtr GetMacroSiteTypeName (Int4 macro_site_type);
NLM_EXTERN ValNodePtr GetSiteTypeList (void);



NLM_EXTERN FieldTypePtr GetFromFieldFromFieldPair (FieldPairTypePtr fieldpair);
NLM_EXTERN FieldTypePtr GetToFieldFromFieldPair (FieldPairTypePtr fieldpair);
NLM_EXTERN FieldPairTypePtr BuildFieldPairFromFromField (FieldTypePtr field_from);
NLM_EXTERN Uint1 FieldTypeFromAECRAction (AECRActionPtr action);
NLM_EXTERN Uint1 FieldTypeChoiceFromFieldPairTypeChoice (Uint1 field_pair_choice);
NLM_EXTERN Int2 FeatureTypeFromFieldType (FieldTypePtr field);
NLM_EXTERN Int4 GetFeatureTypeForRnaType (Int4 rnatype);
NLM_EXTERN int CompareFieldTypes (FieldTypePtr vnp1, FieldTypePtr vnp2);
NLM_EXTERN Boolean AreAECRActionFieldsEqual (AECRActionPtr action1, AECRActionPtr action2);
NLM_EXTERN ValNodePtr GetFieldTypeListFromAECRAction (AECRActionPtr action);
NLM_EXTERN Uint1 GetBiomolForRnaType (Int4 rnatype);
NLM_EXTERN CharPtr GetBiomolNameForRnaType (Int4 rnatype);
NLM_EXTERN void AddAllRNASubtypesToChoiceList (ValNodePtr PNTR field_list);

/* source qual functions */
NLM_EXTERN CharPtr GetSourceQualFromBioSource (BioSourcePtr biop, SourceQualChoicePtr scp, StringConstraintPtr constraint);
NLM_EXTERN ValNodePtr GetMultipleSourceQualsFromBioSource (BioSourcePtr biop, SourceQualChoicePtr scp, StringConstraintPtr constraint);
NLM_EXTERN CharPtr GetQualFromFeature (SeqFeatPtr sfp, FeatureFieldPtr field, StringConstraintPtr scp);
NLM_EXTERN CharPtr GetQualFromFeatureEx (SeqFeatPtr sfp, FeatureFieldPtr field, StringConstraintPtr scp, BatchExtraPtr batch_extra);
NLM_EXTERN Boolean SetQualOnFeature (SeqFeatPtr sfp, FeatureFieldPtr field, StringConstraintPtr scp, CharPtr value, Uint2 existing_text);
NLM_EXTERN Boolean RemoveQualFromFeature (SeqFeatPtr sfp, FeatureFieldPtr field, StringConstraintPtr scp);
NLM_EXTERN Boolean SetSourceQualInBioSource (BioSourcePtr biop, SourceQualChoicePtr scp, StringConstraintPtr constraint, CharPtr value, Uint2 existing_text);
NLM_EXTERN Boolean RemoveSourceQualFromBioSource (BioSourcePtr biop, SourceQualChoicePtr scp, StringConstraintPtr constraint);
NLM_EXTERN Boolean AllowSourceQualMulti (SourceQualChoicePtr s);

NLM_EXTERN ValNodePtr SourceQualValsFromBioSourcePtr (BioSourcePtr biop);
NLM_EXTERN BioSourcePtr BioSourceFromSourceQualVals (ValNodePtr fields);
NLM_EXTERN CharPtr GetDBxrefFromBioSource (BioSourcePtr biop, CharPtr db_name);
NLM_EXTERN Boolean SetDBxrefForBioSource (BioSourcePtr biop, CharPtr db_name, CharPtr str, Uint2 existing_text);
NLM_EXTERN Boolean RemoveDBxrefForBioSource (BioSourcePtr biop, CharPtr db_name, StringConstraintPtr scp);

/* RNA-field functions */
NLM_EXTERN ValNodePtr GetRNATypeList (void);
NLM_EXTERN ValNodePtr GetRnaFieldList (void);
NLM_EXTERN CharPtr GetNameForRnaField (Int4 rnafield);
NLM_EXTERN CharPtr SummarizeRnaType (RnaFeatTypePtr rt);
NLM_EXTERN FeatureFieldPtr FeatureFieldFromRnaQual (RnaQualPtr rq);
NLM_EXTERN RnaQualPtr RnaQualFromFeatureField (FeatureFieldPtr ffp);
NLM_EXTERN Boolean SetRNARefProductString (RnaRefPtr rrp, StringConstraintPtr scp, CharPtr new_val, Uint2 existing_text);
NLM_EXTERN CharPtr GetRNARefProductString (RnaRefPtr rrp, StringConstraintPtr scp);
NLM_EXTERN CharPtr GetRNAProductString (SeqFeatPtr sfp, StringConstraintPtr scp);
NLM_EXTERN Boolean SetRNAProductString (SeqFeatPtr sfp, StringConstraintPtr scp, CharPtr new_val, Uint2 existing_text);
NLM_EXTERN Boolean RemoveRNAProductString (SeqFeatPtr sfp, StringConstraintPtr scp);
NLM_EXTERN Boolean SettmRNATagPeptide (RnaRefPtr rrp, StringConstraintPtr scp, CharPtr new_val, Uint2 existing_text);
NLM_EXTERN CharPtr GettmRNATagPeptide (RnaRefPtr rrp, StringConstraintPtr scp);
NLM_EXTERN Boolean SetncRNAClass (RnaRefPtr rrp, StringConstraintPtr scp, CharPtr new_val, Uint2 existing_text);
NLM_EXTERN CharPtr GetncRNAClass (RnaRefPtr rrp, StringConstraintPtr scp);



NLM_EXTERN CharPtr GetPubFieldLabel (Int4 pub_field);
NLM_EXTERN ValNodePtr GetPubFieldList (void);
NLM_EXTERN CharPtr GetPubFieldFromPub (PubPtr the_pub, Int4 field, StringConstraintPtr scp);
NLM_EXTERN Int4 GetPubMLStatus (PubPtr the_pub);

/* generic string functions */
NLM_EXTERN Boolean SetStringValue (CharPtr PNTR existing_val, CharPtr new_val, Uint2 existing_text);
NLM_EXTERN Boolean RemoveValNodeStringMatch (ValNodePtr PNTR list, StringConstraintPtr scp);
NLM_EXTERN Boolean SetStringsInValNodeStringList (ValNodePtr PNTR list, StringConstraintPtr scp, CharPtr new_val, Uint2 existing_text);


NLM_EXTERN Boolean IsStringConstraintEmpty (StringConstraintPtr scp);
NLM_EXTERN Boolean DoesSingleStringMatchConstraint (CharPtr str, StringConstraintPtr scp);
NLM_EXTERN Boolean DoesStringMatchConstraint (CharPtr str, StringConstraintPtr scp);
NLM_EXTERN Boolean RemoveStringConstraintPortionFromString (CharPtr PNTR str, StringConstraintPtr scp);
NLM_EXTERN Boolean IsSourceConstraintEmpty (SourceConstraintPtr scp);
NLM_EXTERN Boolean DoesBiosourceMatchConstraint (BioSourcePtr biop, SourceConstraintPtr scp);
NLM_EXTERN Boolean IsSequenceConstraintEmpty (SequenceConstraintPtr constraint);
NLM_EXTERN Boolean IsPublicationConstraintEmpty (PublicationConstraintPtr constraint);
NLM_EXTERN Boolean IsFieldConstraintEmpty (FieldConstraintPtr constraint);
NLM_EXTERN Boolean IsCDSGeneProtQualConstraintEmpty (CDSGeneProtQualConstraintPtr constraint);
NLM_EXTERN Boolean IsLocationConstraintEmpty (LocationConstraintPtr lcp);
NLM_EXTERN Boolean DoesObjectMatchConstraintChoiceSet (Uint1 choice, Pointer data, ConstraintChoiceSetPtr csp);
NLM_EXTERN Boolean DoesSeqIDListMeetStringConstraint (SeqIdPtr sip, StringConstraintPtr string_constraint);
NLM_EXTERN ValNodePtr FreeObjectList (ValNodePtr vnp);
NLM_EXTERN ValNodePtr GetObjectListForAECRAction (SeqEntryPtr sep, AECRActionPtr action);
NLM_EXTERN ValNodePtr GetObjectListForAECRActionEx (SeqEntryPtr sep, AECRActionPtr action, BatchExtraPtr batch_extra);
NLM_EXTERN ValNodePtr GetObjectListForFieldType (Uint1 field_type, SeqEntryPtr sep);
NLM_EXTERN Int4 DoApplyActionToObjectList (ApplyActionPtr action, ValNodePtr object_list, Boolean also_change_mrna, StringConstraintPtr scp);
NLM_EXTERN Int4 DoApplyActionToObjectListEx (ApplyActionPtr action, ValNodePtr object_list, Boolean also_change_mrna, StringConstraintPtr scp, BatchExtraPtr batch_extra);
NLM_EXTERN Int4 DoEditActionToObjectList (EditActionPtr action, ValNodePtr object_list, Boolean also_change_mrna);
NLM_EXTERN Int4 DoEditActionToObjectListEx (EditActionPtr action, ValNodePtr object_list, Boolean also_change_mrna, BatchExtraPtr batch_extra);
NLM_EXTERN Int4 DoConvertActionToObjectList (ConvertActionPtr action, ValNodePtr object_list, Boolean also_change_mrna, StringConstraintPtr scp);
NLM_EXTERN Int4 DoConvertActionToObjectListEx (ConvertActionPtr action, ValNodePtr object_list, Boolean also_change_mrna, StringConstraintPtr scp, BatchExtraPtr batch_extra);
NLM_EXTERN Int4 DoCopyActionToObjectList (CopyActionPtr action, ValNodePtr object_list, Boolean also_change_mrna, StringConstraintPtr scp);
NLM_EXTERN Int4 DoCopyActionToObjectListEx (CopyActionPtr action, ValNodePtr object_list, Boolean also_change_mrna, StringConstraintPtr scp, BatchExtraPtr batch_extra);
NLM_EXTERN Int4 DoSwapActionToObjectList (SwapActionPtr action, ValNodePtr object_list, Boolean also_change_mrna, StringConstraintPtr scp);
NLM_EXTERN Int4 DoRemoveActionToObjectList (RemoveActionPtr action, ValNodePtr object_list, Boolean also_change_mrna, StringConstraintPtr scp);
NLM_EXTERN Int4 DoParseActionToObjectList (AECRParseActionPtr action, ValNodePtr object_list, Boolean also_change_mrna, StringConstraintPtr scp);
NLM_EXTERN Int4 DoParseActionToObjectListEx (AECRParseActionPtr action, ValNodePtr object_list, Boolean also_change_mrna, StringConstraintPtr scp, BatchExtraPtr batch_extra);
NLM_EXTERN StringConstraintPtr FindStringConstraintInConstraintSetForField (FieldTypePtr field, ConstraintChoiceSetPtr csp);
NLM_EXTERN StringConstraintPtr FindStringConstraintInConstraintSetForFieldPair (FieldPairTypePtr fieldpair, ConstraintChoiceSetPtr csp);
NLM_EXTERN StringConstraintPtr StringConstraintFromFieldEdit (FieldEditPtr edit);

NLM_EXTERN int LIBCALLBACK SortVnpByObject (VoidPtr ptr1, VoidPtr ptr2);

NLM_EXTERN Boolean IsConversionSupported (Uint2 featdef_from, Uint2 featdef_to);

NLM_EXTERN CharPtr GetTextPortionFromString (CharPtr str, TextPortionPtr text_portion);
NLM_EXTERN Boolean RemoveTextPortionFromString (CharPtr str, TextPortionPtr text_portion);
NLM_EXTERN Boolean IsTextMarkerEmpty (TextMarkerPtr marker);
NLM_EXTERN TextMarkerPtr MakeTextTextMarker (CharPtr text);

NLM_EXTERN Uint2 GetEntityIdFromObject (Uint1 choice, Pointer data);

typedef struct aecrsample {
  FieldTypePtr field;
  CharPtr first_value;
  Int4    num_found;
  Boolean all_same;
} AECRSampleData, PNTR AECRSamplePtr;

NLM_EXTERN AECRSamplePtr AECRSampleFree (AECRSamplePtr sample);

NLM_EXTERN ValNodePtr AECRSampleListFree (ValNodePtr list);
NLM_EXTERN ValNodePtr GetAECRSampleListForSeqEntry (Uint1 field_type, SeqEntryPtr sep);
NLM_EXTERN ValNodePtr GetAECRSampleList (AECRActionPtr act, SeqEntryPtr sep);
NLM_EXTERN AECRSamplePtr GetFieldSampleFromList (ValNodePtr list, FieldTypePtr field);
NLM_EXTERN void GetAECRExistingTextList (Uint1 field_type, SeqEntryPtr sep, FILE *fp);
NLM_EXTERN AECRSamplePtr GetAECRSampleFromObjectList (ValNodePtr object_list, FieldTypePtr field);
NLM_EXTERN AECRSamplePtr GetAECRSampleFromObjectListEx (ValNodePtr object_list, FieldTypePtr field, BatchExtraPtr batch_extra);

NLM_EXTERN AECRSamplePtr GetExistingTextForParseAction (ParseActionPtr action, SeqEntryPtr sep);

NLM_EXTERN void ExportFieldTable (Uint1 field_type, ValNodePtr src_field_list, SeqEntryPtr sep, FILE *fp);

NLM_EXTERN ValNodePtr GetFieldListForFieldType (Uint1 field_type, SeqEntryPtr sep);

NLM_EXTERN ValNodePtr GetSourceQualFieldListFromBioSource (BioSourcePtr biop);
NLM_EXTERN void SortUniqueFieldTypeList (ValNodePtr PNTR field_list);
NLM_EXTERN ValNodePtr LIBCALLBACK FieldTypeListFree (ValNodePtr list);
NLM_EXTERN ValNodePtr LIBCALLBACK FieldTypeListCopy (ValNodePtr orig);

NLM_EXTERN void ApplyMolinfoBlockToSeqEntry (SeqEntryPtr sep, MolinfoBlockPtr mib);

NLM_EXTERN CharPtr GetDescriptorNameFromDescriptorType (Int4 descriptortype);
NLM_EXTERN void AddAllDescriptorsToChoiceList (ValNodePtr PNTR descriptor_type_list);


NLM_EXTERN void ApplyMacroToSeqEntry (SeqEntryPtr sep, ValNodePtr macro, Int4Ptr pNumFields, Int4Ptr pNumFeat);

NLM_EXTERN SeqFeatPtr ApplyOneFeatureToBioseq (BioseqPtr bsp, Uint1 featdef, SeqLocPtr slp, ValNodePtr fields, ValNodePtr src_fields, Boolean add_mrna);

/* for generating text representations of macro objects */
NLM_EXTERN CharPtr SummarizeSourceQual (ValNodePtr field);
NLM_EXTERN CharPtr FeatureFieldLabel (CharPtr feature_name, ValNodePtr field);
NLM_EXTERN Boolean IsFeatureFieldEmpty (FeatureFieldPtr field);
NLM_EXTERN Boolean IsFieldTypeEmpty (FieldTypePtr field);
NLM_EXTERN CharPtr SummarizeFieldType (ValNodePtr vnp);
NLM_EXTERN Boolean AllowFieldMulti (FieldTypePtr field);

NLM_EXTERN CharPtr GetFieldValueForObject (Uint1 choice, Pointer data, FieldTypePtr field, StringConstraintPtr scp);
NLM_EXTERN CharPtr GetFieldValueForObjectEx (Uint1 choice, Pointer data, FieldTypePtr field, StringConstraintPtr scp, BatchExtraPtr batch_extra);
NLM_EXTERN Boolean SetFieldValueForObject (Uint1 choice, Pointer data, FieldTypePtr field, StringConstraintPtr scp, CharPtr value, Uint2 existing_text);
NLM_EXTERN Boolean SetFieldValueForObjectEx (Uint1 choice, Pointer data, FieldTypePtr field, StringConstraintPtr scp, CharPtr value, Uint2 existing_text, BatchExtraPtr batch_extra);
NLM_EXTERN BioseqPtr GetSequenceForObject (Uint1 choice, Pointer data);
NLM_EXTERN ValNodePtr GetMultipleFieldValuesForObject (Uint1 choice, Pointer data, FieldTypePtr field, StringConstraintPtr scp, BatchExtraPtr batch_extra);

typedef enum {
  eTableMatchFeatureID = 1,
  eTableMatchGeneLocusTag,
  eTableMatchProteinID,
  eTableMatchDbxref,
  eTableMatchNucID,
  eTableMatchBioSource,
  eTableMatchSourceQual
} ETableMatchType;


typedef struct matchtype {
  Uint1   choice;
  Pointer data;
  Uint1   match_location;
} MatchTypeData, PNTR MatchTypePtr;

typedef struct tabcolumnconfig {
  MatchTypePtr match_type;
  FieldTypePtr field;
  Uint2 existing_text;
  Boolean skip_blank;
  Boolean match_mrna;
  ValNodePtr constraint;
} TabColumnConfigData, PNTR TabColumnConfigPtr;

NLM_EXTERN MatchTypePtr MatchTypeNew ();
NLM_EXTERN MatchTypePtr MatchTypeFree (MatchTypePtr match_type);

NLM_EXTERN TabColumnConfigPtr TabColumnConfigNew (void);
NLM_EXTERN TabColumnConfigPtr TabColumnConfigFree (TabColumnConfigPtr t);
NLM_EXTERN TabColumnConfigPtr TabColumnConfigCopy (TabColumnConfigPtr orig);
NLM_EXTERN ValNodePtr TabColumnConfigListFree (ValNodePtr columns);
NLM_EXTERN ValNodePtr TabColumnConfigListCopy (ValNodePtr orig);
NLM_EXTERN ValNodePtr ValidateTabTableValues (ValNodePtr table, ValNodePtr columns);
NLM_EXTERN ValNodePtr ValidateFeatureFieldColumnNames (ValNodePtr header_line, ValNodePtr PNTR perr_list);
NLM_EXTERN ValNodePtr FreeObjectTableForTabTable (ValNodePtr table);
NLM_EXTERN ValNodePtr GetObjectTableForTabTable (SeqEntryPtr sep, ValNodePtr table, ValNodePtr columns, ValNodePtr PNTR p_err_list);
NLM_EXTERN ValNodePtr GetSequenceListsForMatchTypeInTabTable (SeqEntryPtr sep, ValNodePtr table, Int4 col, MatchTypePtr match_type, ValNodePtr PNTR p_err_list);
NLM_EXTERN ValNodePtr FreeSequenceLists (ValNodePtr lists);
NLM_EXTERN ValNodePtr ApplyTableValuesToObjectTable (SeqEntryPtr sep, ValNodePtr table, ValNodePtr columns, ValNodePtr obj_table);
NLM_EXTERN ValNodePtr CheckObjTableForRowsThatApplyToTheSameDestination (ValNodePtr obj_table);
NLM_EXTERN ValNodePtr CheckObjTableForExistingText (SeqEntryPtr sep, ValNodePtr table, ValNodePtr columns, ValNodePtr obj_table);

NLM_EXTERN ValNodePtr ApplyTableToFeatures (SeqEntryPtr sep, ValNodePtr table, ValNodePtr columns);
NLM_EXTERN ValNodePtr CheckTableForExistingText (SeqEntryPtr sep, ValNodePtr table, ValNodePtr columns);

NLM_EXTERN SeqFeatPtr GetmRNAForFeature (SeqFeatPtr sfp);
NLM_EXTERN Boolean AdjustmRNAProductToMatchProteinProduct (SeqFeatPtr sfp);
NLM_EXTERN Boolean IsFieldTypeCDSProduct (FieldTypePtr ft);

NLM_EXTERN ValNodePtr GetPublicationTitlesInSep (SeqEntryPtr sep);
NLM_EXTERN ValNodePtr GetBankitCommentsInSep (SeqEntryPtr sep);
NLM_EXTERN ValNodePtr GetPublicationTitlesOnSep (SeqEntryPtr sep);
NLM_EXTERN ValNodePtr GetBankitCommentsOnSep (SeqEntryPtr sep);

NLM_EXTERN BioseqPtr GetRepresentativeBioseqFromBioseqSet (BioseqSetPtr bssp);

NLM_EXTERN ValNodePtr ValNodeCopyPtr (ValNodePtr orig);
NLM_EXTERN SeqLocPtr ParseSimpleSeqLoc (CharPtr str, BioseqPtr bsp);

NLM_EXTERN void FixCapitalizationInString (CharPtr PNTR pTitle, Uint2 capitalization, ValNodePtr   org_names);

NLM_EXTERN Boolean GBBlockIsCompletelyEmpty (GBBlockPtr gb);

NLM_EXTERN CharPtr GetObjectIdString (ObjectIdPtr oip);
NLM_EXTERN Boolean SetObjectIdString (ObjectIdPtr oip, CharPtr value, Uint2 existing_text);

#ifdef __cplusplus 
} 
#endif

#endif
