FROM ubuntu:trusty

RUN apt-get update && apt-get -y install python python-pip python-dev build-essential
RUN pip install numpy pandas
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US:en
ENV LC_ALL en_US.UTF-8

COPY LMSDFDownload6Dec16FinalAll.sdf.converted /main/LMSDFDownload6Dec16FinalAll.sdf.converted
COPY entry.py /main/entry.py

WORKDIR /main/

ENTRYPOINT ["python", "./entry.py", "-d", "/main/LMSDFDownload6Dec16FinalAll.sdf.converted", "-i", "/input/", "-o", "/output/"]


