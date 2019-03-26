Bootstrap: docker
From: buisciii/centos7_base_image:latest

%files
    ./scif_app_recipes/* /opt/
%post
    echo "Install basic development tools"
    yum -y groupinstall "Development Tools"
    yum -y update && yum -y install wget curl openssl-devel geos-devel udunits2-devel libxml2-devel cairo-devel libgit2-devel


    echo "Install python2.7 setuptools and pip"
    yum -y install python-setuptools
    easy_install pip

    echo "Installing SCI-F"
    pip install scif

    echo "Installing FastQC app" && \
    scif install /opt/fastqc_v0.11.7_centos7.scif
    echo "Installing seqtk app" && \
    scif install /opt/seqtk_v1.3_centos7.scif
    echo "Installing NGSQCToolkit app" && \
    scif install /opt/
    echo "Installing trimmomatic app" && \
    scif install /opt/trimmomatic_v0.38_centos7.scif && \
    echo "Installing samtools app" && \
    scif install /opt/samtools_v1.2_centos7.scif && \
    echo "Installing bamutil app" && \
    scif install /opt/bamutil_v1.0.13_centos7.scif
    echo "Installing picard app" && \
    scif install /opt/picard_v1.140_centos7.scif && \
    echo "Installing IGV app" && \
    scif install /opt/IGV_v2.4.19_centos7.scif
    echo "Installing bwa app" && \
    scif install /opt/bwa_v0.7.17_centos7.scif && \
    echo "Installing vep app" && \
    scif install /opt/vep_v95.2_centos7.scif
    echo "Installing GATK app" && \
    scif install /opt/gatk_v3.8_centos7.scif
    echo "Installing vcftools app" && \
    scif install /opt/vcftools_v0.1.12b_centos7.scif && \
    echo "Installing NGSQCToolkit app" && \
    scif install /opt/NGSQCToolkit_v2.3.3_centos7.scif && \
    echo "Installing seq_crumbs app" && \
    scif install /opt/seqcrumbs_v0.1.9_centos7.scif && \
    echo "Installing spades app" && \
    scif install /opt/spades_v3.8.0_centos7.scif && \
    echo "Installing quast app" && \
    scif install /opt/quast_v5.0.0_centos7.scif && \
    echo "Installing srst2 app" && \
    scif install /opt/srst2_v0.1.8_centos7.scif && \
    echo "Installing multiqc app" && \
    scif install /opt/multiqc_v1.4_centos7.scif

    # Executables must be exported for nextflow, if you use their singularity native integration.
    # It would be cool to use $SCIF_APPBIN_bwa variable, but it must be set after PATH variable, because I tried to use it here and in %environment without success.
    find /scif/apps -maxdepth 2 -name "bin" | while read in; do echo "export PATH=\${PATH}:$in" >> $SINGULARITY_ENVIRONMENT;done

    find /scif/apps -maxdepth 2 -name "lib" | while read in; do echo "export LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:$in" >> $SINGULARITY_ENVIRONMENT ;done

    if [[ ":$PATH:" == *":/scif/apps/snppipeline:"* ]];then

        export CLASSPATH=/scif/apps/varscan/varscan-2.3.9/varscan-2.3.9.jar:$CLASSPATH >> $SINGULARITY_ENVIRONMENT
        export CLASSPATH=/scif/apps/picard/picard.jar:$CLASSPATH >> $SINGULARITY_ENVIRONMENT
        export CLASSPATH=/scif/apps/gatk/gatk-3.8/GenomeAnalysisTK.jar:$CLASSPATH >> $SINGULARITY_ENVIRONMENT

    fi

%runscript
    exec scif "$@"
