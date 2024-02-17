#!/bin/bash

mv *.sum ./sum/
mv *.vex ./vex/
mv *.key ./key/

rm *.oms *.ke *.wa *.hb *.cd *.tv2d *.flag *.ho *.yg
rm *.vex2
rm fort.*

echo "Moved important files and deleted annoying files"

