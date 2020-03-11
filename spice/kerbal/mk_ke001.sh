rm -v ke001.bsp
mkspk -setup moho_mkspk.txt    -input moho.txt   -output ke001.bsp
mkspk -setup eve_mkspk.txt     -input eve.txt    -output ke001.bsp -append
mkspk -setup gilly_mkspk.txt   -input gilly.txt  -output ke001.bsp -append
mkspk -setup kerbin_mkspk.txt  -input kerbin.txt -output ke001.bsp -append

