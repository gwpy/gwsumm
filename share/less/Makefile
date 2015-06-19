all: gwsumm-l1.css gwsumm-h1.css gwsumm-v1.css gwsumm-c1.css gwsumm-all.css

gwsumm-l1.css: gwsumm-l1.less gwsumm.less
	lessc --global-var="ifo=l1" gwsumm.less -x > ../css/$@;

gwsumm-h1.css: gwsumm-h1.less gwsumm.less
	lessc --global-var="ifo=h1" gwsumm.less -x > ../css/$@;

gwsumm-v1.css: gwsumm-v1.less gwsumm.less
	lessc --global-var="ifo=v1" gwsumm.less -x > ../css/$@;

gwsumm-c1.css: gwsumm-c1.less gwsumm.less
	lessc --global-var="ifo=c1" gwsumm.less -x > ../css/$@;

gwsumm-all.css: gwsumm-all.less gwsumm.less
	lessc --global-var="ifo=all" gwsumm.less -x > ../css/$@;

distribute:
	cp gwsumm-l1.css ~/public_html/html/ligo/gwsumm.css
	gsiscp gwsumm-h1.css ldas-pcdev2.ligo-wa.caltech.edu:~/public_html/html/ligo/gwsumm.css
	gsiscp gwsumm-all.css ldas-pcdev2.ligo.caltech.edu:~/public_html/html/ligo/gwsumm.css

distribute-detchar:
	gsiscp gwsumm-l1.css ~detchar/public_html/html/ligo/gwsumm.css
	gsiscp gwsumm-h1.css ldas-pcdev2.ligo-wa.caltech.edu:~detchar/public_html/html/ligo/gwsumm.css
	gsiscp gwsumm-all.css ldas-pcdev2.ligo.caltech.edu:~detchar/public_html/html/ligo/gwsumm.css
