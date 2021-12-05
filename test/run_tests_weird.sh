MATSIZE=16384

for RESTART in {1..50}
do
    echo ./build/krylov.exe -m import -n $MATSIZE -i $MATSIZE -r $(expr $RESTART \* 5) -k 4 -o 2 >> ./results/n${MATSIZE}/results_2_k4.txt
    ./build/krylov.exe -m import -n $MATSIZE -i $MATSIZE -r $(expr $RESTART \* 5) -k 4 -o 2 >> ./results/n${MATSIZE}/results_2_k4.txt
    echo "" >> ./results/n${MATSIZE}/results_2_k4.txt
done


echo ./build/krylov.exe -m import -n $MATSIZE -i $MATSIZE -r 4 -k 4 -o 4 >> ./results/n${MATSIZE}/results_4_k4.txt
./build/krylov.exe -m import -n $MATSIZE -i $MATSIZE -r 4 -k 4 -o 4 >> ./results/n${MATSIZE}/results_4_k4.txt
echo "" >> ./results/n${MATSIZE}/results_4_k4.txt

for RESTART in {1..50}
do
    echo ./build/krylov.exe -m import -n $MATSIZE -i $MATSIZE -r $(expr $RESTART \* 10) -k 4 -o 4 >> ./results/n${MATSIZE}/results_4_k4.txt
    ./build/krylov.exe -m import -n $MATSIZE -i $MATSIZE -r $(expr $RESTART \* 10) -k 4 -o 4 >> ./results/n${MATSIZE}/results_4_k4.txt
    echo "" >> ./results/n${MATSIZE}/results_4_k4.txt
done