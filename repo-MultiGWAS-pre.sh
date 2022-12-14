dirm=$1/multiGWAS-pre
echo "Copy multiGWAS precompiledin " $dirm
read

# root
mkdir -p $dirm
cp README.md $dirm
cp INSTALL.sh $dirm
cp UNINSTALL.sh $dirm
cp multiGWAS-gui.png $dir m

cp -a main $dirm
cp -a examples $dirm

# install
mkdir $dirm/install
cp -a install/repo $dirm/install
cp install/install-* $dirm/install

# opt
opt=$dirm/opt
mkdir $opt
cp -a opt/tools $opt/
cp -a opt/Rlibs $opt/
cp -a opt/jre   $opt/

# docs/supplements
sups=docs/supplements
mkdir -p $dirm/$sups
cp $sups/supp*.pdf $dirm/$sups
cp $sups/supp*.html $dirm/$sups

