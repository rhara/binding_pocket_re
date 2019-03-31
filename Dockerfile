FROM deepchemio/deepchem:latest

RUN apt-get update
RUN apt-get upgrade -y
RUN apt-get install vim git bzip2 -y

WORKDIR /root

