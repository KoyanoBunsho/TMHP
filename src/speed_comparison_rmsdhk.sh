#!/bin/sh

g++ specific_rmsdhk.cpp -o specific_rmsdhk -std=c++14 -lstdc++fs -Wall -Wextra -O3 -mtune=native -march=native -mfpmath=both -Werror 
g++ specific_fast_rmsdhk.cpp -o specific_fast_rmsdhk -std=c++14 -lstdc++fs -Wall -Wextra -O3 -mtune=native -march=native -mfpmath=both -Werror
g++ specific_fast_rmsdhk_postpro.cpp -o specific_fast_rmsdhk_postpro -std=c++14 -lstdc++fs -Wall -Wextra -O3 -mtune=native -march=native -mfpmath=both -Werror

g++ specific_rmsdhk_more_data.cpp -o specific_rmsdhk_more_data -std=c++14 -lstdc++fs -Wall -Wextra -O3 -mtune=native -march=native -mfpmath=both -Werror 
g++ specific_fast_rmsdhk_more_data.cpp -o specific_fast_rmsdhk_more_data -std=c++14 -lstdc++fs -Wall -Wextra -O3 -mtune=native -march=native -mfpmath=both -Werror
g++ specific_fast_rmsdhk_more_data_postpro.cpp -o specific_fast_rmsdhk_more_data_postpro -std=c++14 -lstdc++fs -Wall -Wextra -O3 -mtune=native -march=native -mfpmath=both -Werror

for k in {2..5}
do
echo "-----$k-----"
echo "-----RMSDhk-----"
./specific_rmsdhk $k
echo "-----Linear Time RMSDhk-----"
./specific_fast_rmsdhk $k
./specific_fast_rmsdhk_postpro $k

echo "-----More data-----"

echo "-----RMSDhk-----"
./specific_rmsdhk_more_data $k
echo "-----Linear Time RMSDhk-----"
./specific_fast_rmsdhk_more_data $k
./specific_fast_rmsdhk_more_data_postpro $k
done
