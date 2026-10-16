#!/usr/bin/env Rscript
# End-to-end synthetic synthesis tests, including deliberately corrupted outputs.
# Run from repository root: Rscript tests/public_release/test-public-synthesis.R
library(data.table)
library(yaml)
root <- tempfile("public-synthesis-"); dir.create(root)
write_table <- function(x, path) { dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE); fwrite(x, path) }
run <- function(script, cfg, expect_success = TRUE, failure_pattern = NULL) {
  log <- tempfile(tmpdir = root)
  status <- system2("Rscript", c(shQuote(script), "--config", shQuote(cfg)), stdout = log, stderr = log)
  if ((status == 0) != expect_success) stop(paste(readLines(log), collapse = "\n"))
  if (!is.null(failure_pattern) && !any(grepl(failure_pattern, readLines(log), fixed = TRUE))) stop("Unexpected validation failure: ", paste(readLines(log), collapse = "\n"))
  cat(if (expect_success) "PASS" else "PASS rejected corruption", script, "\n")
}
mag <- file.path(root, "mag"); dyn <- file.path(root, "dynamics")
wmag <- function(x, path) write_table(x, file.path(mag, path))
wdyn <- function(x, path) write_table(x, file.path(dyn, path))
for (p in c(file.path(mag, "00_quality_control/10_audit_candidate_annotations.complete"), file.path(dyn, "00_quality_control/04_make_figures_and_report.complete"))) { dir.create(dirname(p), recursive = TRUE); writeLines("synthetic complete", p) }
metrics <- c("gut_mags_represented_by_members", "unique_parsed_mags", "unique_member_proteins", "unique_representatives", "cross_mag_representative_clusters", "representatives_with_annotation_row", "representatives_with_functional_annotation", "representatives_without_annotation_row")
wmag(data.table(metric = metrics, value = c(6, 6, 100, 40, 20, 30, 25, 10), percent_of_relevant_denominator = c(100,100,100,100,50,75,62.5,25)), "00_quality_control/Table_02_full_catalogue_qc_summary.csv")
wmag(data.table(size_group = c("1", "2-5", "6-20", "21-100"), clusters = c(10,20,30,40), unknown_or_unannotated_fraction = c(.1,.4,.8,.6)), "02_unknown_by_cluster_size/Table_04_unknown_annotation_by_cluster_size.csv")
wmag(data.table(odds_ratio = .9, lower_95 = .5, upper_95 = 1.2, p_value = .2, status = "fitted"), "02_unknown_by_cluster_size/singleton_unknown_enrichment_model.csv")
wmag(data.table(ontology = c("KO", "EC"), term = c("K04028", "1.2.1.87"), omnibus_p_adj_bh = c(.8,.9)), "04_longitudinal_functional_potential/Table_06_longitudinal_functional_models.csv.gz")
tr <- CJ(candidate_id = c("KO:K04028", "EC:1.2.1.87"), Subject = paste0("synthetic_subject_", 1:3), Timepoint = c("T1", "T2", "T6", "T10"))
tr[, sampleID := paste(Subject, Timepoint)]; tr[, potential := seq_len(.N)/100]
wmag(tr, "06_targeted_candidate_validation/candidate_participant_trajectories.csv.gz")
con <- data.table(MAG = paste0("synthetic_MAG_",1:6), candidate_1_carried = c(T,T,T,T,T,F), candidate_2_carried = c(F,F,T,T,T,T))
con[, shared_carrier := candidate_1_carried & candidate_2_carried]
wmag(con, "10_candidate_annotation_audit/Table_24_candidate_carrier_MAG_concordance.csv")
wmag(data.table(representative_id = c("synthetic_rep", "synthetic_rep"), representative_mag = c("synthetic_MAG_3", "synthetic_MAG_3"), carries_both_candidate_terms_same_representative = TRUE), "10_candidate_annotation_audit/Table_23_candidate_representative_protein_clusters.csv")
wmag(data.table(candidate_id = c("KO:K04028", "EC:1.2.1.87"), representative_clusters_carrying_both_terms = 1, carrier_MAGs_with_both_terms_on_same_representative = 1, carrier_union_MAGs = 6, shared_carrier_MAGs = 3), "10_candidate_annotation_audit/Table_25_candidate_curation_summary.csv")
wmag(data.table(display_class = c("KO:K04028 only", "Shared by both candidates", "EC:1.2.1.87 only"), fill_class = c("K04028 only", "Shared", "EC 1.2.1.87 only"), carrier_mags = c(2,3,1)), "reproducibility/Figure_11_plot_data.csv")
wmag(data.table(component = "shared_carriers", model_family = "absolute_component_temporal_pattern", carrier_mags = 3, p_adj_bh_within_family = .8), "08_shared_unique_carrier_decomposition/Table_17_component_temporal_models.csv")
wmag(data.table(contrast_id = "T2_vs_T6", before_timepoint = "T2", after_timepoint = "T6", ratio_after_vs_before = .7, p_adj_bh_across_45_timepoint_contrasts = .9), "09_shared_component_followup/Table_18_shared_component_timepoint_contrasts.csv")
wmag(data.table(contrast_id = "T2_vs_T6", genus = paste0("g__Synthetic", 1:10), mean_delta = -(10:1)/100, absolute_change_share = (10:1)/55), "09_shared_component_followup/Table_21_shared_component_genus_drivers.csv.gz")
wdyn(data.table(layer = c("taxonomic_MAG", "KO"), metric = "hellinger", permanova_r_squared = c(.1,.05), permanova_p_value = c(.4,.8), permanova_p_adj_bh_functional_family = c(NA,.9)), "02_multivariate_analysis/Table_27_global_PERMANOVA_and_dispersion_results.csv")
t <- data.table(layer = "KO", subject_id = rep(paste0("synthetic_subject_",1:3), each = 3), taxonomic_standardized_distance = seq(.1,.9,.1), functional_standardized_distance = c(.9,.6,.8,.4,.7,.3,.5,.2,.1))
wdyn(t, "reproducibility/Figure_18_plot_data.csv")
wdyn(data.table(layer = "KO", model_family = "taxonomic_functional_transition_coupling"), "03_longitudinal_stability/Table_30_matched_taxonomic_functional_stability_models.csv")
wdyn(data.table(layer = "KO", sensitivity_family = c("participant_spearman_transition_coupling", "participant_median_relative_movement"), participant_summary = c(-.2,.7), positive_direction_participants = c(1,1), participants = 3, p_adj_bh_within_sensitivity_family = .8, ratio = c(NA,.7)), "03_longitudinal_stability/Table_31_participant_level_transition_sensitivity.csv")
wdyn(data.table(layer = "KO", subject_id = paste0("synthetic_subject_",1:3), recovery_fraction_from_peak = c(.2,-.3,.4)), "reproducibility/Figure_19_plot_data.csv")
cfg <- read_yaml("config/manuscript_functional_synthesis.yml"); cfg$inputs <- list(mag_functional_results_dir = mag, functional_dynamics_results_dir = dyn); cfg$outputs$output_dir <- file.path(root, "functional"); cfg$outputs$png_dpi <- 60; cfg$outputs$table_preview_dpi <- 60
cfg_path <- file.path(root,"functional.yml"); write_yaml(cfg,cfg_path)
builder <- "analysis/manuscript_functional_synthesis/00_build_main_figure_and_tables_public.R"
validator <- "analysis/manuscript_functional_synthesis/tests/validate_outputs_public.R"
run(builder,cfg_path); run(validator,cfg_path)
table2 <- fread(file.path(cfg$outputs$output_dir,"main/Main_Table_02_key_statistical_results.csv"))
stopifnot(!any(grepl("Roseburia|higher at T6|machine precision", unlist(table2))), any(grepl("Synthetic1",table2$Analysis)))
# Changed inputs must be accepted; a copied, mismatching result must be rejected.
panel_path <- file.path(cfg$outputs$output_dir,"reproducibility/Main_Figure_01_panel_c.csv"); original <- fread(panel_path)
bad <- copy(original); bad$carrier_mags[1] <- bad$carrier_mags[1]+1; fwrite(bad,panel_path); run(validator,cfg_path,FALSE,"Carrier counts disagree"); fwrite(original,panel_path)
bad <- copy(original); bad$segment_label[1] <- "wrong"; fwrite(bad,panel_path); run(validator,cfg_path,FALSE,"Carrier diagram labels"); fwrite(original,panel_path)
# A source-evidence mismatch must stop synthesis; zero overlap is also legitimate.
curation_path <- file.path(mag,"10_candidate_annotation_audit/Table_25_candidate_curation_summary.csv")
curation <- fread(curation_path); bad_curation <- copy(curation)
bad_curation$representative_clusters_carrying_both_terms <- 2L; fwrite(bad_curation,curation_path)
cfg$outputs$output_dir <- file.path(root,"functional-bad-curation"); write_yaml(cfg,cfg_path)
run(builder,cfg_path,FALSE,"Curation summary disagrees")
curation$representative_clusters_carrying_both_terms <- 0L
curation$carrier_MAGs_with_both_terms_on_same_representative <- 0L; fwrite(curation,curation_path)
wmag(data.table(representative_id=character(), representative_mag=character(), carries_both_candidate_terms_same_representative=logical()),"10_candidate_annotation_audit/Table_23_candidate_representative_protein_clusters.csv")
cfg$outputs$output_dir <- file.path(root,"functional-zero-overlap"); write_yaml(cfg,cfg_path)
run(builder,cfg_path); run(validator,cfg_path)
# NetMoss: synthetic taxa, counts, scores and edge identities unlike the study.
net <- file.path(root,"netmoss"); n <- 11; taxa <- paste0("g__Synthetic",1:n)
write_table(data.table(taxon_names = taxa, NetMoss_Score = seq(.05,.95,length.out=n), p.val = seq(.01,.8,length.out=n), p.adj = c(.01,.02,rep(.8,n-2))), file.path(net,"tables/netmoss_scores.csv"))
source_manifest <- list(status="COMPLETED_EXPLORATORY", group_a=list(timepoints=c("T1","T2","Ta"), independent_profiles=3,raw_samples=6), group_b=list(timepoints=c("T8","T9","T10"),independent_profiles=4,raw_samples=9),retained_taxa=n,scenario="synthetic", generated_at="synthetic",git_commit="synthetic")
dir.create(file.path(net,"reproducibility")); write_yaml(source_manifest,file.path(net,"reproducibility/netmoss_manifest.yml"))
for (group in c("a","b")) {
 abundance <- matrix(seq_len(n*(if(group=="a") 3 else 4))/1000,nrow=n,dimnames=list(taxa,paste0("synthetic_profile_",seq_len(if(group=="a") 3 else 4))))
 write_table(data.frame(taxon=taxa,abundance,check.names=FALSE),file.path(net,paste0("prepared/group_",group,"_abundance.csv")))
 mat <- diag(n); dimnames(mat)<-list(taxa,taxa)
 mat[1,2]<-mat[2,1]<-if(group=="a") .6 else -.7
 mat[2,3]<-mat[3,2]<-.3
 write_table(data.frame(taxon=taxa,mat,check.names=FALSE),file.path(net,paste0("networks/group_",group,"_network.csv")))
}
gut_path<-file.path(root,"gut.csv"); write_table(data.frame(sampleID="synthetic_sample",setNames(as.list(rep(.01,n)),paste0("synthetic_MAG_",1:n)),check.names=FALSE),gut_path)
tax_path<-file.path(root,"taxonomy.csv"); write_table(data.table(Genome=paste0("synthetic_MAG_",1:n),genus=taxa,phylum=paste0("p__Synthetic",rep(1:5,length.out=n))),tax_path)
ncfg <- read_yaml("config/netmoss_manuscript_synthesis.yml"); ncfg$inputs <- list(netmoss_results_dir=net,gut_abundance_file=gut_path,taxonomy_file=tax_path); ncfg$outputs$output_dir<-file.path(root,"net-synthesis"); ncfg$figure$dpi<-60
ncfg_path<-file.path(root,"net.yml");write_yaml(ncfg,ncfg_path)
nbuilder<-"analysis/netmoss_manuscript_synthesis/00_build_main_figure_and_tables_public.R"; nvalidator<-"analysis/netmoss_manuscript_synthesis/tests/validate_outputs_public.R"
run(nbuilder,ncfg_path); run(nvalidator,ncfg_path)
caption<-paste(readLines(file.path(ncfg$outputs$output_dir,"main/Main_Figure_01_caption.txt")),collapse=" ")
stopifnot(grepl("11 scored genera",caption,fixed=TRUE),grepl("3 and 4 independent profiles",caption,fixed=TRUE))
edges_path<-file.path(ncfg$outputs$output_dir,"reproducibility/Main_Figure_01_panel_d_edges.csv"); original_edges<-fread(edges_path)
bad<-copy(original_edges); bad$taxon_1[1]<-"g__Wrong"; fwrite(bad,edges_path);run(nvalidator,ncfg_path,FALSE,"Displayed edge identities disagree");fwrite(original_edges,edges_path)
bad<-copy(original_edges);bad$retained[1]<-!bad$retained[1];fwrite(bad,edges_path);run(nvalidator,ncfg_path,FALSE,"Edge values or retention flags disagree");fwrite(original_edges,edges_path)
# Zero edges, fewer phyla and no FDR findings must still render and validate.
for (group in c("a", "b")) {
 mat <- diag(n); dimnames(mat) <- list(taxa, taxa)
 write_table(data.frame(taxon = taxa, mat, check.names = FALSE), file.path(net, paste0("networks/group_", group, "_network.csv")))
}
scores_path <- file.path(net, "tables/netmoss_scores.csv")
scores <- fread(scores_path); scores$p.adj <- .8; fwrite(scores, scores_path)
small_tax <- fread(tax_path); small_tax$phylum <- "p__Synthetic"; fwrite(small_tax, tax_path)
ncfg$outputs$output_dir <- file.path(root,"net-empty"); write_yaml(ncfg,ncfg_path)
run(nbuilder,ncfg_path); run(nvalidator,ncfg_path)
cat("Synthetic synthesis checks passed.\n")
