#! /usr/bin/env bats

@test "TESTS: rnftools et2roc" {
	./command_line/rnftools_et2roc.sh
}
@test "TESTS: rnftools es2et" {
	./command_line/rnftools_es2et.sh
}
@test "TESTS: rnftools curesim2rnf" {
	./command_line/rnftools_curesim2rnf.sh
}
@test "TESTS: rnftools wgsim2rnf" {
	./command_line/rnftools_wgsim2rnf.sh
}
@test "TESTS: rnftools sam2roc" {
	./command_line/rnftools_sam2roc.sh
}
@test "TESTS: rnftools check" {
	./command_line/rnftools_check.sh
}
@test "TESTS: rnftools art2rnf" {
	./command_line/rnftools_art2rnf.sh
}
@test "TESTS: rnftools publication" {
	./command_line/rnftools_publication.sh
}
@test "TESTS: rnftools sam2es" {
	./command_line/rnftools_sam2es.sh
}
@test "TESTS: rnftools merge" {
	./command_line/rnftools_merge.sh
}
@test "TESTS: rnftools mason2rnf" {
	./command_line/rnftools_mason2rnf.sh
}
@test "TESTS: rnftools dwgsim2rnf" {
	./command_line/rnftools_dwgsim2rnf.sh
}
@test "TESTS: rnftools liftover" {
	./command_line/rnftools_liftover.sh
}

#################################################

@test "TESTS: snakemake: all simulators SE" {
	./snakemake/01_*/run.sh -p
}

@test "TESTS: snakemake: all simulators PE" {
	./snakemake/02_*/run.sh -p
}

@test "TESTS: snakemake: many files SE" {
	./snakemake/03_*/run.sh -p
}
	
@test "TESTS: snakemake: many files PE" {
	./snakemake/04_*/run.sh -p
}

@test "TESTS: snakemake: zero reads SE" {
	./snakemake/05_*/run.sh -p
}

@test "TESTS: snakemake: zero reads PE" {
	./snakemake/06_*/run.sh -p
}

#################################################
#################################################

@test "EXAMPLES: 01_tutorial/01" {
	../examples/01_tutorial/01_*/run.sh -p --cores
}

#################################################

@test "EXAMPLES: 01_tutorial/02/01" {
	../examples/01_tutorial/02_*/01_*/run.sh -p --cores
}

@test "EXAMPLES: 01_tutorial/02/01 RNF" {
	./_test_rnf.sh ../examples/01_tutorial/02_*/01_*/
}

###

@test "EXAMPLES: 01_tutorial/02/02" {
	../examples/01_tutorial/02_*/02_*/run.sh -p --cores
}

@test "EXAMPLES: 01_tutorial/02/02 RNF" {
	./_test_rnf.sh ../examples/01_tutorial/02_*/02_*/
}

###

@test "EXAMPLES: 01_tutorial/02/03" {
	../examples/01_tutorial/02_*/03_*/run.sh -p --cores
}

@test "EXAMPLES: 01_tutorial/02/03 RNF" {
	./_test_rnf.sh ../examples/01_tutorial/02_*/03_*/
}

###

@test "EXAMPLES: 01_tutorial/02/04" {
	../examples/01_tutorial/02_*/04_*/run.sh -p --cores
}

@test "EXAMPLES: 01_tutorial/02/04 RNF" {
	./_test_rnf.sh ../examples/01_tutorial/02_*/04_*/
}

###

@test "EXAMPLES: 01_tutorial/02/05" {
	../examples/01_tutorial/02_*/05_*/run.sh -p --cores
}

@test "EXAMPLES: 01_tutorial/02/05 RNF" {
	./_test_rnf.sh ../examples/01_tutorial/02_*/05_*/
}

###

@test "EXAMPLES: 01_tutorial/02/06" {
	../examples/01_tutorial/02_*/06_*/run.sh -p --cores
}

@test "EXAMPLES: 01_tutorial/02/06 RNF" {
	./_test_rnf.sh ../examples/01_tutorial/02_*/06_*/
}


#################################################

@test "EXAMPLES: 01_tutorial/03/01" {
	../examples/01_tutorial/03_*/01_*/run.sh -p --cores
}

@test "EXAMPLES: 01_tutorial/03/01 RNF" {
	./_test_rnf.sh ../examples/01_tutorial/03_*/01_*/bams/
}

###

@test "EXAMPLES: 01_tutorial/03/02" {
	../examples/01_tutorial/03_*/02_*/run.sh -p --cores
}

@test "EXAMPLES: 01_tutorial/03/02 RNF" {
	./_test_rnf.sh ../examples/01_tutorial/03_*/02_*/bams/
}

###

@test "EXAMPLES: 01_tutorial/03/03" {
	../examples/01_tutorial/03_*/03_*/run.sh -p --cores
}

@test "EXAMPLES: 01_tutorial/03/03 RNF" {
	./_test_rnf.sh ../examples/01_tutorial/03_*/03_*/bams/
}


#################################################

@test "EXAMPLES: 02" {
	../examples/02_*/run.sh -p --cores
}

@test "EXAMPLES: 02 RNF" {
	./_test_rnf.sh ../examples/02_*/
}

###

@test "EXAMPLES: 04" {
	../examples/04_*/run.sh -p --cores
}

#################################################

@test "DOCUMENTATION" {
	cd ../docs
	make clean
	mkdir -p _static
	sphinx-build -W -b html -d _build/doctrees . _build/html
}

#################################################


