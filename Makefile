run:
	nvidia-docker run -t -i -v ~/data:/data -v $$PWD:/supp -p 8888:8888 deepchemio/deepchem bash
