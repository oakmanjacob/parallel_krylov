MATSIZE=16384

RESTART=75

for KVAL in {1..10}
do
echo ./build/krylov.exe -m import -n $MATSIZE -i $MATSIZE -r $RESTART -k $KVAL -o 2 >> ./results/n${MATSIZE}/results_2_r${RESTART}_t1.txt
./build/krylov.exe -m import -n $MATSIZE -i $MATSIZE -r $RESTART -k $KVAL -o 2 >> ./results/n${MATSIZE}/results_2_r${RESTART}_t1.txt
echo "" >> ./results/n${MATSIZE}/results_2_r${RESTART}_t1.txt
done

for TVAL in {2..8}
do
    for KVAL in {1..10}
    do
    echo ./build/krylov.exe -m import -n $MATSIZE -i $MATSIZE -r $RESTART -k $KVAL -o 3 -t $TVAL -p 10 >> ./results/n${MATSIZE}/results_3_r${RESTART}_t${TVAL}.txt
    ./build/krylov.exe -m import -n $MATSIZE -i $MATSIZE -r $RESTART -k $KVAL -o 3 -t $TVAL -p 10 >> ./results/n${MATSIZE}/results_3_r${RESTART}_t${TVAL}.txt
    echo "" >> ./results/n${MATSIZE}/results_3_r${RESTART}_t${TVAL}.txt
    done
done