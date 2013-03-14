sphinx-apidoc -o . ../cobra
rm *oven*rst
rm *db_tools*rst
rm *omics*rst
rm *internal*rst
DOCDIR=../cobra/documentation
rm $DOCDIR/doctrees/*
rm $DOCDIR/html/*html
rm $DOCDIR/html/*js
rm $DOCDIR/html/_static/*
rm $DOCDIR/html/_sources/*
DOCDIR=$DOCDIR/html/_modules/
rm $DOCDIR/*html
rm $DOCDIR/*/*html
rm $DOCDIR/*/*/*html
