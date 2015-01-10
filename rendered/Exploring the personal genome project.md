
# Exploring the Personal Genome Project

This notebook demonstrates working with sample genomics data stored as publicly
accessible BigQuery datasets.

### In this notebook you will

* Explore the Personal Genome Project genomics datasets available in BigQuery
* Use the `%%bq_sql` statement to write and execute SQL statements within the
notebook
* Extract data from BigQuery and create a local dataset that can be manipulated
in Python
* Refine and pivot your local dataset via the Pandas Python library
* Visualize different aspects of your dataset via the Matplotlib Python library


### Attribution
This notebook is based on the [Google
Genomics](https://cloud.google.com/genomics/) BigQuery examples here:
* [Personal Genome Project](https://github.com/googlegenomics/bigquery-
examples/tree/master/pgp)

----


    import gcp.bigquery as bq

# Exploring the Personal Genome Project dataset

First let's look at the table of phenotypes that we want to explore here. We can
use the `gcp.bigquery` Python package to fetch the table schema. We can see a
number of fields that provide details for a given participant's genome, such as
the existence of certain diseases/conditions (e.g., has_Asthma) as well as the
participant's heritage (e.g., Maternal_grandfather_Country_of_origin).


    phenotypes = bq.table('google.com:biggene:pgp.phenotypes')
    phenotypes.schema()




<table><tr><th>name</th><th>data_type</th><th>mode</th><th>description</th></tr><tr><td>Participant</td><td>STRING</td><td>NULLABLE</td><td></td></tr><tr><td>Year_of_birth</td><td>STRING</td><td>NULLABLE</td><td></td></tr><tr><td>Which_statement_best_describes_you</td><td>STRING</td><td>NULLABLE</td><td></td></tr><tr><td>Severe_disease_or_rare_genetic_trait</td><td>STRING</td><td>NULLABLE</td><td></td></tr><tr><td>Disease_trait_Onset</td><td>STRING</td><td>NULLABLE</td><td></td></tr><tr><td>Disease_trait_Rarity</td><td>STRING</td><td>NULLABLE</td><td></td></tr><tr><td>Disease_trait_Severity</td><td>STRING</td><td>NULLABLE</td><td></td></tr><tr><td>Disease_trait_Relative_enrollment</td><td>STRING</td><td>NULLABLE</td><td></td></tr><tr><td>Disease_trait_Diagnosis</td><td>STRING</td><td>NULLABLE</td><td></td></tr><tr><td>Disease_trait_Genetic_confirmation</td><td>STRING</td><td>NULLABLE</td><td></td></tr><tr><td>Disease_trait_Documentation</td><td>STRING</td><td>NULLABLE</td><td></td></tr><tr><td>Sex_Gender</td><td>STRING</td><td>NULLABLE</td><td></td></tr><tr><td>Race_ethnicity</td><td>STRING</td><td>NULLABLE</td><td></td></tr><tr><td>Maternal_grandmother_Country_of_origin</td><td>STRING</td><td>NULLABLE</td><td></td></tr><tr><td>Paternal_grandmother_Country_of_origin</td><td>STRING</td><td>NULLABLE</td><td></td></tr><tr><td>Paternal_grandfather_Country_of_origin</td><td>STRING</td><td>NULLABLE</td><td></td></tr><tr><td>Maternal_grandfather_Country_of_origin</td><td>STRING</td><td>NULLABLE</td><td></td></tr><tr><td>Enrollment_of_relatives</td><td>STRING</td><td>NULLABLE</td><td></td></tr><tr><td>Enrollment_of_older_individuals</td><td>STRING</td><td>NULLABLE</td><td></td></tr><tr><td>Enrollment_of_parents</td><td>STRING</td><td>NULLABLE</td><td></td></tr><tr><td>Enrolled_relatives_Monozygotic_Identical_twins</td><td>STRING</td><td>NULLABLE</td><td></td></tr><tr><td>Enrolled_relatives_Parents</td><td>STRING</td><td>NULLABLE</td><td></td></tr><tr><td>Enrolled_relatives_Siblings_Fraternal_twins</td><td>STRING</td><td>NULLABLE</td><td></td></tr><tr><td>Enrolled_relatives_Children</td><td>STRING</td><td>NULLABLE</td><td></td></tr><tr><td>Enrolled_relatives_Grandparents</td><td>STRING</td><td>NULLABLE</td><td></td></tr><tr><td>Enrolled_relatives_Grandchildren</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>Enrolled_relatives_Aunts_Uncles</td><td>STRING</td><td>NULLABLE</td><td></td></tr><tr><td>Enrolled_relatives_Nephews_Nieces</td><td>STRING</td><td>NULLABLE</td><td></td></tr><tr><td>Enrolled_relatives_Half_siblings</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>Enrolled_relatives_Cousins_or_more_distant</td><td>STRING</td><td>NULLABLE</td><td></td></tr><tr><td>Enrolled_relatives_Not_genetically_related_e_g_husband_wife</td><td>STRING</td><td>NULLABLE</td><td></td></tr><tr><td>Are_all_your_enrolled_relatives_linked_to_your_PGP_profile</td><td>STRING</td><td>NULLABLE</td><td></td></tr><tr><td>Have_you_uploaded_genetic_data_to_your_PGP_participant_profile</td><td>STRING</td><td>NULLABLE</td><td></td></tr><tr><td>Have_you_used_the_PGP_web_interface_to_record_a_designated_proxy</td><td>STRING</td><td>NULLABLE</td><td></td></tr><tr><td>Have_you_uploaded_health_record_data_using_our_Google_Health_or_Microsoft_Healthvault_interfaces</td><td>STRING</td><td>NULLABLE</td><td></td></tr><tr><td>Uploaded_health_records_Update_status</td><td>STRING</td><td>NULLABLE</td><td></td></tr><tr><td>Uploaded_health_records_Extensiveness</td><td>INTEGER</td><td>NULLABLE</td><td></td></tr><tr><td>Blood_sample</td><td>STRING</td><td>NULLABLE</td><td></td></tr><tr><td>Saliva_sample</td><td>STRING</td><td>NULLABLE</td><td></td></tr><tr><td>Microbiome_samples</td><td>STRING</td><td>NULLABLE</td><td></td></tr><tr><td>Tissue_samples_from_surgery</td><td>STRING</td><td>NULLABLE</td><td></td></tr><tr><td>Tissue_samples_from_autopsy</td><td>STRING</td><td>NULLABLE</td><td></td></tr><tr><td>Month_of_birth</td><td>STRING</td><td>NULLABLE</td><td></td></tr><tr><td>Anatomical_sex_at_birth</td><td>STRING</td><td>NULLABLE</td><td></td></tr><tr><td>Maternal_grandmother_Race_ethnicity</td><td>STRING</td><td>NULLABLE</td><td></td></tr><tr><td>Maternal_grandfather_Race_ethnicity</td><td>STRING</td><td>NULLABLE</td><td></td></tr><tr><td>Paternal_grandmother_Race_ethnicity</td><td>STRING</td><td>NULLABLE</td><td></td></tr><tr><td>Paternal_grandfather_Race_ethnicity</td><td>STRING</td><td>NULLABLE</td><td></td></tr><tr><td>has_NA</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Sciatica</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Allergic_contact_dermatitis</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Urinary_tract_infection_UTI</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Endometriosis</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Ovarian_cysts</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Impacted_tooth</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Dental_cavities</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Gallstones</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Chronic_tonsillitis</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Allergic_rhinitis</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Asthma</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Hyperopia_Farsightedness</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Iron_deficiency_anemia</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Migraine_with_aura</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Lactose_intolerance</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Cervical_cancer</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Osteoarthritis</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Achilles_tendonitis</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Dandruff</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Keloids</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Skin_tags</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Gingivitis</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Temporomandibular_joint_TMJ_disorder</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Canker_sores_oral_ulcers</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Gastroesophageal_reflux_disease_GERD</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Hiatal_hernia</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Nonalcoholic_fatty_liver_disease_NAFLD</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Hypertension</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Hemorrhoids</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Meniere_s_disease</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Tinnitus</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_High_triglycerides_hypertriglyceridemia</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Thyroid_nodule_s</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Hashimoto_s_thyroiditis</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Colon_polyps</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Scoliosis</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Lichen_planus</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Cafe_au_lait_spots</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Myopia_Nearsightedness</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Floaters</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Hypothyroidism</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_High_cholesterol_hypercholesterolemia</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Presbyopia</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Trigger_finger</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Dupuytren_s_contracture</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Plantar_fasciitis</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Fibromyalgia</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Eczema</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Fibrocystic_breast_disease</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Premature_ventricular_contractions</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Cardiac_arrhythmia</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Restless_legs_syndrome</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Cluster_headaches</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Carpal_tunnel_syndrome</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Other_peripheral_neuropathy</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Uterine_fibroids</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Flatfeet</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Inguinal_hernia</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Age_related_hearing_loss</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Non_melanoma_skin_cancer</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Hair_loss_includes_female_and_male_pattern_baldness</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Acne</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Rosacea</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Geographic_tongue</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Irritable_bowel_syndrome_IBS</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Spinal_stenosis</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Kidney_stones</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Male_infertility</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Diverticulosis</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Deviated_septum</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Diabetes_mellitus</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_type_2</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Astigmatism</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Peptic_ulcer_stomach_or_duodenum</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Appendicitis</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Age_related_cataract</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Bone_spurs</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Rheumatoid_arthritis</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Hemochromatosis</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Graves_disease</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Congenital_heart_defect</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Chronic_bronchitis</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Mitral_valve_prolapse</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Epilepsy</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Hemophilia</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Melanoma</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Chronic_liver_disease_and_cirrhosis</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Age_related_macular_degeneration</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Sensorineural_hearing_loss_or_congenital_deafness</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Hyperhidrosis_excessive_sweating</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Fissured_tongue</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Nasal_polyps</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Chronic_sinusitis</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Bunions</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Osteoporosis</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Strabismus</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Amyotrophic_lateral_sclerosis_ALS</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Frozen_shoulder</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Single_transverse_palmar_crease_simian_crease</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Dermatographia</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Rotator_cuff_tear</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Lipoma</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Polycystic_ovary_syndrome_PCOS</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Hereditary_motor_and_sensory_neuropathy_includes_Charcot_Marie_Tooth_disease_and_HNPP</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Varicose_veins</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Retinal_detachment</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Prostate_cancer</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Spermatocele</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Varicocele</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Essential_tremor</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Hereditary_thrombophilia_includes_Factor_V_Leiden_and_Prothrombin_G20210A</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Breast_cancer</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Celiac_disease</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Multiple_sclerosis_MS</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Migraine_without_aura</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Ehlers_Danlos_syndrome</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Bartholin_s_cyst</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Tennis_elbow</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Angina</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Atrial_fibrillation</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Idiopathic_thrombocytopenic_purpura_ITP</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Heart_block</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Pilonidal_cyst</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Psoriasis</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Bundle_branch_block</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Chondromalacia_patella_CMP</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Rectal_prolapse</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Raynaud_s_phenomenon</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Dry_eye_syndrome</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Other_thrombophilia_includes_antiphospholipid_syndrome</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Tongue_tie_ankyloglossia</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Marfan_syndrome</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Color_blindness</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Breast_fibroadenoma</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Crohn_s_disease</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Thyroid_cancer</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Non_Hodgkin_lymphoma</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Female_infertility</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Sjogren_s_syndrome_Sicca_syndrome</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Benign_prostatic_hypertrophy_BPH</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Osgood_Schlatter_disease</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Postural_kyphosis</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Parkinson_s_disease</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Narcolepsy</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Pernicious_anemia_a_k_a_Addison_Biermer_anemia</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Trigeminal_neuralgia</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Non_Celiac_Gluten_Sensitivity</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Gout</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Deep_vein_thrombosis_DVT</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Hereditary_hemorrhagic_telangiectasia_also_known_as_Osler_Weber_Rendu_syndrome</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Ulcerative_colitis</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Polydactyly</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Uterine_prolapse</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Gilbert_syndrome</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Spina_bifida</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Long_QT_Syndrome</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_type_1</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Myocardial_infarction_heart_attack</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Congestive_heart_failure</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Glaucoma</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Chronic_tension_headaches_15_days_per_month</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_at_least_6_months</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Barrett_s_esophagus</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Peyronie_s_disease</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Diabetic_retinopathy</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Lung_cancer</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Acute_liver_failure</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Sick_sinus_syndrome_includes_tachy_brady_syndrome</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Traumatic_cataract</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Hypospadias</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Stomach_cancer</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Colon_cancer</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Rectal_cancer</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Pancreatic_cancer</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Endometrial_cancer</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Ovarian_cancer</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Bladder_cancer</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Kidney_cancer</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Leukemia</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Polycythemia_vera</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Essential_thrombocythemia</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Neurofibromatosis</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Brain_cancer</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Wolff_Parkinson_White_WPW_Syndrome</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Cleft_palate</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Recurrent_sleep_paralysis</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Folate_deficiency_anemia</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Syndactyly_webbing_of_digits</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Alopecia_areata</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Scheuermann_s_kyphosis</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Bell_s_palsy</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Growth_hormone_deficiency</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Pulmonary_embolism</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Other_cardiomyopathy_including_ARVD</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Von_Willebrand_disease</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Polycystic_kidney_disease</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Stroke</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Acute_kidney_failure</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Otosclerosis</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Lupus</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Arnold_Chiari_malformation</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Congenital_nystagmus</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Emphysema</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Chronic_Obstructive_Pulmonary_Disease_COPD</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Congenital_clubfoot_equinovarus</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Alpha_1_antitrypsin_deficiency</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Sickle_cell_trait_carrier</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Infantile</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_juvenile</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_and_presenile_cataract</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Huntington_s_disease</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Porphyria</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Hidradenitis_suppurativa</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Aortic_aneurysm</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Dilated_cardiomyopathy</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Urethral_diverticulum</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Diverticulitis_Urinary_Tract</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Hypertrophic_cardiomyopathy</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Keratoconus</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Central_serous_retinopathy</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Developmental_dysplasia_of_the_hip</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Retinitis_pigmentosa</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Cleft_uvula</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Muscular_dystrophy</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Congenital_ichthyosis</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Autoimmune_hemolytic_anemia</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Infantile_pyloric_stenosis</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Hypertensive_retinopathy</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_G6PD_deficiency</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Other_aneurysm</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr><tr><td>has_Congenital_hydronephrosis</td><td>BOOLEAN</td><td>NULLABLE</td><td></td></tr></table>



Let's look at a few rows from the table for a few phenotype attributes of
interest.


    phenotypes.sample(fields=['Participant',
                              'Sex_Gender',
                              'Year_of_birth',
                              'Maternal_grandmother_Country_of_origin',
                              'Maternal_grandfather_Country_of_origin'])




<div>Number of rows: 5</div><div>Query job ID  : job_Z3GjYLseTkx8vw8z2GCxoSwFgmY</div>
<table><tr><th>Year_of_birth</th><th>Participant</th><th>Maternal_grandmother_Country_of_origin</th><th>Sex_Gender</th><th>Maternal_grandfather_Country_of_origin</th></tr><tr><td>1973</td><td>huBFBBD8</td><td>United States</td><td>Male</td><td>United States</td></tr><tr><td>1986</td><td>hu5880D9</td><td>United States</td><td>Female</td><td>United States</td></tr><tr><td>1987</td><td>hu9367D1</td><td>Poland</td><td>Male</td><td>Austria</td></tr><tr><td>1960</td><td>hu57850F</td><td>United States</td><td>Female</td><td>United States</td></tr><tr><td>&nbsp;</td><td>hu1069B0</td><td>China</td><td>Female</td><td>China</td></tr></table>



## Querying

Let's dig further into the phenotypes table and get some statistics on the
occurrence of asthma among the participants present within the dataset. The
`%%bq_sql` statement allows us to write SQL within our notebook and execute it
within BigQuery.


    %%bq_sql asthma
    SELECT Participant,
           IFNULL(Sex_Gender, 'Unknown') AS gender,
           IF(has_Asthma IS NULL, 0, 1) AS has_asthma
    FROM $phenotypes

Our query is now stored in a variable called `asthma`. We can execute the query
(by sending to BigQuery) by calling `asthma.results()`.  Additionally, let's
convert the results to a [Pandas dataframe](http://pandas.pydata.org/pandas-
docs/dev/generated/pandas.DataFrame.html) and display just the first 5 rows from
the results.


    asthma.results().to_dataframe()[:5]




<div style="max-height:1000px;max-width:1500px;overflow:auto;">
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Participant</th>
      <th>gender</th>
      <th>has_asthma</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td> huBFBBD8</td>
      <td>   Male</td>
      <td> 0</td>
    </tr>
    <tr>
      <th>1</th>
      <td> hu5880D9</td>
      <td> Female</td>
      <td> 0</td>
    </tr>
    <tr>
      <th>2</th>
      <td> hu9367D1</td>
      <td>   Male</td>
      <td> 0</td>
    </tr>
    <tr>
      <th>3</th>
      <td> hu57850F</td>
      <td> Female</td>
      <td> 0</td>
    </tr>
    <tr>
      <th>4</th>
      <td> hu1069B0</td>
      <td> Female</td>
      <td> 0</td>
    </tr>
  </tbody>
</table>
</div>



### Composing Queries

Now that we've seen the structure of the rows generated by our `$asthma` query,
let's summarize the results by computing the percentage of the overall dataset
having asthma. We can accomplish this by nesting our previous query within an
outer `SELECT` statement. The `%%bq_sql` statement allows us to embed our
previous `$asthma` query like so:


    %%bq_sql
    SELECT ROUND(AVG(has_asthma) * 100, 2) AS percent_has_asthma FROM $asthma




<div>Number of rows: 1</div><div>Query job ID  : job_OVxGHy7VlePTpXqrqreFD9CaUqc</div>
<table><tr><th>percent_has_asthma</th></tr><tr><td>9.1</td></tr></table>



Notice the reference to the previous query we defined via the `$asthma`
placeholder. The resulting query is a nested SQL statement.  We can see the
full, nested query that was generated:


    print asthma.sql

    SELECT Participant,
           IFNULL(Sex_Gender, 'Unknown') AS gender,
           IF(has_Asthma IS NULL, 0, 1) AS has_asthma
    FROM [google.com:biggene:pgp.phenotypes]


This is handy because we can build up our large, complex query one piece at a
time and validate that the intermediate results we're getting along the way.

Let's further refine our `$asthma` query and compute the average occurrence
broken down by gender this time.


    %%bq_sql
    SELECT gender,
           ROUND(AVG(has_asthma) * 100, 2) AS percent_has_asthma
    FROM $asthma
    GROUP BY gender




<div>Number of rows: 4</div><div>Query job ID  : job_Iq_SS48lqcB8-4cSFh3eO9S6FJQ</div>
<table><tr><th>gender</th><th>percent_has_asthma</th></tr><tr><td>Male</td><td>6.8</td></tr><tr><td>Female</td><td>12.03</td></tr><tr><td>Transmasculine </td><td>0.0</td></tr><tr><td>Unknown</td><td>13.92</td></tr></table>



Looks like there are roughly twice as many female (versus male) genomes with the
asthma phenotype in this dataset.

### Python Pandas Integration

Since our `$asthma` query returns a result set that is small enough to easily
fit into memory (at ~2k rows), let's create a local dataset from our query by
populating a [Pandas dataframe](http://pandas.pydata.org/pandas-
docs/dev/generated/pandas.DataFrame.html) object and saving it to the `df`
variable.


    df = asthma.results().to_dataframe()

This dataframe allows us to slice and dice our data in a number of ways and all
of the computation will happen locally. Let's try getting the number of genomes
per gender via our dataframe.


    df.groupby('gender').sum()




<div style="max-height:1000px;max-width:1500px;overflow:auto;">
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>has_asthma</th>
    </tr>
    <tr>
      <th>gender</th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>Female</th>
      <td> 103</td>
    </tr>
    <tr>
      <th>Male</th>
      <td>  85</td>
    </tr>
    <tr>
      <th>Transmasculine </th>
      <td>   0</td>
    </tr>
    <tr>
      <th>Unknown</th>
      <td>  11</td>
    </tr>
  </tbody>
</table>
</div>



So that shows us the total number of participants annotated to have asthma,
broken down by gender. Let's see if we can reproduce the percent_has_asthma
column that we computed above via BigQuery.


    df.groupby('gender').mean() * 100




<div style="max-height:1000px;max-width:1500px;overflow:auto;">
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>has_asthma</th>
    </tr>
    <tr>
      <th>gender</th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>Female</th>
      <td> 12.032710</td>
    </tr>
    <tr>
      <th>Male</th>
      <td>  6.800000</td>
    </tr>
    <tr>
      <th>Transmasculine </th>
      <td>  0.000000</td>
    </tr>
    <tr>
      <th>Unknown</th>
      <td> 13.924051</td>
    </tr>
  </tbody>
</table>
</div>



Great! So for datasets that easily fit into memory, we can easily aggregate our
results via the Pandas dataframe using an approach very similar to our SQL
statement.

## Visualization

Every Pandas dataframe can be used for visualization by calling the
[`dataframe.plot()`](http://pandas.pydata.org/pandas-
docs/dev/visualization.html) method.  Let's plot the counts of asthma occurrence
broken down by gender.


    tally = df.groupby('gender').sum()
    tally.plot(kind='pie', y='has_asthma', legend=False)




    <matplotlib.axes.AxesSubplot at 0x109522210>




![png](output_30_1.png)



    tally.plot(kind='bar', y='has_asthma', legend=False)




    <matplotlib.axes.AxesSubplot at 0x109897a90>




![png](output_31_1.png)

