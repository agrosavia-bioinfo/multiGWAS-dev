# Create users MultiGWAS directory for upload to github copying dirs/files from MultiGWAS-dev

OUTDIR=$1
MGDIR=$OUTDIR/MultiGWAS

if [ -d "$MGDIR" ]; then
	mv $MGDIR $OUTDIR/old-MultiGWAS
fi

find examples -name gwas.log|xargs rm -f
find examples -name gwas.errors|xargs rm -f
find examples -name out-full-tetra-Filters-ADD|xargs rm -rf

mkdir $MGDIR

cp INSTALL_Rlibs.R INSTALL.sh README.md UNINSTALL.sh $MGDIR
cp -r sources docs tools examples  $MGDIR

mkdir $MGDIR/opt/
cp opt/* $MGDIR/opt
mkdir /$MGDIR/opt/Rlibs



