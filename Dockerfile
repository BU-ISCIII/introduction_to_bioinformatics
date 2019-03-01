FROM centos:latest

COPY ./scif_app_recipes/* /opt/


RUN echo "Install basic development tools" && \
    yum -y groupinstall "Development Tools" && \
    yum -y update && yum -y install wget curl && \
    echo "Install python2.7 setuptools and pip" && \
    yum -y install python-setuptools && \
    easy_install pip && \
    echo "Installing SCI-F" && \
    pip install scif ipython

RUN echo "Installing FastQC app" && \
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
    echo "Installing multiqc app" && \
    scif install /opt/multiqc_v1.4_centos7.scif

#ENTRYPOINT ["/opt/docker-entrypoint.sh"]
#CMD ["scif"]

RUN find /scif/apps -maxdepth 2 -name "bin" | while read in; do echo "export PATH=\$PATH:$in" >> /etc/bashrc;done 
RUN if [ -z "${LD_LIBRARY_PATH-}" ]; then echo "export LD_LIBRARY_PATH=/usr/local/lib" >> /etc/bashrc;fi
RUN find /scif/apps -maxdepth 2 -name "lib" | while read in; do echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:$in" >> /etc/bashrc;done
