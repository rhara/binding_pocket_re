build:
	nvidia-docker build -t rhara/deepchem:0.1 .

run:
	nvidia-docker run -t -i -v ~/data:/data -v $$PWD:/supp -p 8888:8888 rhara/deepchem:0.1 bash
