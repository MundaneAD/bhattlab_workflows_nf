
process binning_prep {
	tag "BINNING_PREP on $sample_id"

	input:
	tuple val(sample_id), path(reads), path(contigs)

	output:
	tuple val(sample_id), path(contigs), path("${sample_id}.depth.txt"), path("align_${sample_id}.bam"), emit: bin_prep
	path "versions.yml", emit: versions

	script:
	
	"""
	# make an index for the contigs and align the reads to get the depth
	mkdir -p ./idx_${sample_id}
	mv $contigs ./idx_${sample_id}
	bwa index ./idx_${sample_id}/${contigs}
	
	bwa mem -t $task.cpus ./idx_${sample_id}/${contigs} ${reads[0]} ${reads[1]} \
		| samtools sort --threads $task.cpus > align_${sample_id}.bam
	samtools sort -@ $task.cpus align_${sample_id}.bam

	# what about singles? guess we ignore them for now?
	# seems that like jgi_summarize can take several bam files?

	jgi_summarize_bam_contig_depths --outputDepth ${sample_id}.depth.txt \
		--pairedContigs ${sample_id}.paired.txt --minContigLength 1000 \
		--minContigDepth 1 --percentIdentity 50 ./align_${sample_id}.bam

	mv ./idx_${sample_id}/${contigs} ./
	# add versions here
	bwa 2> bwa.version || true
	jgi_summarize_bam_contig_depths -h 2> jgi.version || true
	cat <<-END_VERSIONS > versions.yml
	"${task.process}":
	    bwa: \$( cat bwa.version | grep Version | sed -e "s/Version: //g" )
	    samtools: \$( samtools --version | head -n 1 | sed -e "s/samtools //g" )
	    metabat2: \$( head -n 1 jgi.version | sed -e "s/jgi_summarize_bam_contig_depths //g" | sed -e "s/ (Bioconda).*//g" )
	END_VERSIONS
	"""
}


process binning_prep_lr_bam {
	tag "BINNING_PREP bam on $sample_id"

	input:
	tuple val(sample_id), path(info)

	output:
	tuple val(sample_id), path("${info[0]}"), path("${info[1]}"), path("${sample_id}.bam"), emit: bin_bam
	path "versions.yml", emit: versions

	script:
	"""
	minimap2 -x map-ont -t ${task.cpus} \
		-a ${info[1]} ${info[0]} | samtools sort --threads ${task.cpus} > ${sample_id}.bam

	cat <<-END_VERSIONS > versions.yml
	"${task.process}":
	    minimap2: \$( minimap2 --version )
	    samtools: \$( samtools --version | head -n 1 | sed -e "s/samtools //g" )
	END_VERSIONS
	"""

}

process binning_prep_lr {
	tag "BINNING_PREP on $sample_id"

	input:
	tuple val(sample_id), path(reads), path(contigs), path(bam)

	output:
	tuple val(sample_id), path(contigs), path("${sample_id}.depth.txt"), path(bam), emit: bin_prep
	path "versions.yml", emit: versions

	script:
	
	"""

	jgi_summarize_bam_contig_depths --outputDepth ${sample_id}.depth.txt \
		--pairedContigs ${sample_id}.paired.txt --minContigLength 1000 \
		--minContigDepth 1 --percentIdentity 50 ${bam}

	jgi_summarize_bam_contig_depths -h 2> jgi.version || true
	cat <<-END_VERSIONS > versions.yml
	"${task.process}":
	    samtools: \$( samtools --version | head -n 1 | sed -e "s/samtools //g" )
	    metabat2: \$( head -n 1 jgi.version | sed -e "s/jgi_summarize_bam_contig_depths //g" | sed -e "s/ (Bioconda).*//g" )
	END_VERSIONS
	"""
}

process binning_prep_hybrid {
    tag "BINNING_PREP_HYBRID for $sample_id"

    input:
    tuple val(sample_id), path(mags), path(short_reads), path(long_reads)

    output:
    tuple val(sample_id), path(mags), path("${sample_id}.depth.txt"), path("${sample_id}_short.bam"), path("${sample_id}_long.bam"), emit: bin_prep
    path "versions.yml", emit: versions

    script:
    """
    # Prepare depth information from the opera-ms .fasta file headers
    mkdir -p ./output_${sample_id}
    mv $mags ./output_${sample_id}

    # Parse contig headers to extract depth information
    awk '/^>/ {
        split(\$0, arr, " ");
        contig_id = substr(arr[1], 2);  # Remove '>' from contig name
        cov_short = substr(arr[3], 11);  # Extract cov_short value
        cov_long = substr(arr[4], 10);  # Extract cov_long value
        print contig_id "\\t" cov_short "\\t" cov_long;
    }' ./output_${sample_id}/${mags} > ${sample_id}.depth.txt

    # Align short reads to contigs using minimap2
    minimap2 -ax sr -t ${task.cpus} ./output_${sample_id}/${mags} ${short_reads[0]} ${short_reads[1]} | \
        samtools sort --threads ${task.cpus} -o ${sample_id}_short.bam

    # Align long reads to contigs using minimap2
    minimap2 -ax map-ont -t ${task.cpus} ./output_${sample_id}/${mags} ${long_reads} | \
        samtools sort --threads ${task.cpus} -o ${sample_id}_long.bam

    # Index BAM files
    samtools index ${sample_id}_short.bam
    samtools index ${sample_id}_long.bam

    # Add versions information
    minimap2 --version > minimap2.version
    samtools --version | head -n 1 > samtools.version
    awk --version | head -n 1 > awk.version

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: $(cat minimap2.version)
        samtools: $(cat samtools.version)
        awk: $(cat awk.version)
    END_VERSIONS
    """
}
