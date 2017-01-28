sed -i "s/UFixed8/N0f8/g" $*
sed -i "s/UFixed10/N6f10/g" $*
sed -i "s/UFixed12/N4f12/g" $*
sed -i "s/UFixed14/N2f14/g" $*
sed -i "s/UFixed16/N0f16/g" $*
sed -i "s/UFixed/Normed/g" $*
# For types from ColorTypes
sed -i "s/U8/N0f8/g" $*
sed -i "s/U16/N0f16/g" $*
