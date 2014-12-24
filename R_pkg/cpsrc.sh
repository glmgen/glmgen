
PREFIX=/usr0/home/vsadhana/code/glmgen/glmgen/
FROM=$PREFIX/R_pkg/glmgen/src/
TO=$PREFIX/c_lib/glmgen/src/

if [ "$1" = "ctor" ]; then
  echo "Copying source from c_lib to R_pkg"
  FROM=$PREFIX/c_lib/glmgen/src/
  TO=$PREFIX/R_pkg/glmgen/src/  
else
  echo "Copying source from R_pkg to c_lib"
fi;

/bin/cp -f $FROM/tf/*.h $FROM/tf/*.c $TO/tf/
/bin/cp -f $FROM/utils/*.h $FROM/utils/*.c $TO/utils/

    



