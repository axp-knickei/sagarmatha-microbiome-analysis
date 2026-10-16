# ============================================
# Makefile for the legacy configured mixed-model pipeline
# ============================================
# Usage:
#   make setup      — Create Conda environment
#   make data       — Prepare raw data
#   make analyze    — Run the configured legacy mixed model
#   make results    — Generate tables & figures
#   make all        — Full pipeline
#   make clean      — Remove outputs

SHELL := /bin/bash
.ONESHELL:

.PHONY: setup data analyze results hpc-candidates alpha-diversity beta-env beta-test beta-diversity all clean

setup:
	conda env create -f environment.yml
	@echo "Environment 'glmm-hpc' created. Activate with: conda activate glmm-hpc"

data:
	conda run -n glmm-hpc Rscript scripts/01_data_prep.R

analyze:
	conda run -n glmm-hpc Rscript scripts/02_run_glmm.R

results:
	conda run -n glmm-hpc Rscript scripts/03_results.R

hpc-candidates:
	conda run -n glmm-hpc Rscript scripts/summarize_hpc_candidate_signals.R

alpha-diversity:
	conda run -n glmm-hpc Rscript scripts/alpha_diversity/01_prepare_alpha_diversity_inputs.R
	conda run -n glmm-hpc Rscript scripts/alpha_diversity/02_calculate_alpha_diversity_vegan.R
	conda run -n glmm-hpc python scripts/alpha_diversity/03_calculate_alpha_diversity_scikitbio.py
	conda run -n glmm-hpc Rscript scripts/alpha_diversity/04_compare_alpha_diversity_tools.R
	conda run -n glmm-hpc Rscript scripts/alpha_diversity/05_plot_alpha_diversity_publication.R

beta-env:
	bash scripts/validate_beta_diversity_environment.sh

beta-test:
	Rscript analysis/beta_diversity/09_create_smoke_test_data.R --config config/beta_diversity_test.yml
	bash scripts/run_beta_diversity_local.sh --config config/beta_diversity_test.yml --force

beta-diversity:
	bash scripts/run_beta_diversity_local.sh --config config/beta_diversity.yml

all: data analyze results

clean:
	rm -f results/tables/*
	rm -f results/figures/*
	rm -f results/logs/*
	rm -f data/processed/*
