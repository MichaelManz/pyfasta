FROM ubuntu:18.04

RUN apt-get update \
  && apt-get install -y python-pip python-dev wget \
  && cd /usr/local/bin \
  && ln -s /usr/bin/python python \
  && pip install --upgrade pip

RUN apt-get install -y libgsl23 \
    && wget https://www.tbi.univie.ac.at/RNA/download/ubuntu/ubuntu_18_04/viennarna_2.4.14-1_amd64.deb \
    && dpkg -i viennarna_2.4.14-1_amd64.deb

VOLUME /data

ADD ./pyfasta /pyfasta/pyfasta
ADD ./setup.* /pyfasta/
ADD ./README.rst /pyfasta/
ADD ./CHANGELOG.txt /pyfasta/
RUN pip install /pyfasta

ADD ./scripts /pyfasta/scripts

CMD ["/bin/bash", "/pyfasta/scripts/processing.sh"]