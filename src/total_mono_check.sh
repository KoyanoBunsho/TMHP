#!/bin/sh

g++ specific_rmsdhk_total_mono_check.cpp -o specific_rmsdhk_total_mono_check -std=c++14 -lstdc++fs -Wall -Wextra -O3 -mtune=native -march=native -mfpmath=both -Werror 
g++ specific_rmsdhk_total_mono_check_more_data.cpp -o specific_rmsdhk_total_mono_check_more_data -std=c++14 -lstdc++fs -Wall -Wextra -O3 -mtune=native -march=native -mfpmath=both -Werror 
for k in {2}
do
echo "-----$k-----"
echo "-----RMSDhk-----"
./specific_rmsdhk_total_mono_check $k
echo "-----More data-----"
./specific_rmsdhk_total_mono_check_more_data $k
done
