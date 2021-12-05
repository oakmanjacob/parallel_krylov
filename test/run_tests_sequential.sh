MATSIZE=1024
for KVAL in {1..10}
do
echo ./build/krylov.exe -m import -n $MATSIZE -i $MATSIZE -r $MATSIZE -k $KVAL -o 1 >> ./results/n${MATSIZE}/results_1.txt
./build/krylov.exe -m import -n $MATSIZE -i $MATSIZE -r $MATSIZE -k $KVAL -o 1 >> ./results/n${MATSIZE}/results_1.txt
echo "" >> ./results/n${MATSIZE}/results_1.txt
done

MATSIZE=4096
for KVAL in {1..10}
do
echo ./build/krylov.exe -m import -n $MATSIZE -i $MATSIZE -r $MATSIZE -k $KVAL -o 1 >> ./results/n${MATSIZE}/results_1.txt
./build/krylov.exe -m import -n $MATSIZE -i $MATSIZE -r $MATSIZE -k $KVAL -o 1 >> ./results/n${MATSIZE}/results_1.txt
echo "" >> ./results/n${MATSIZE}/results_1.txt
done

MATSIZE=16384
for KVAL in {1..10}
do
echo ./build/krylov.exe -m import -n $MATSIZE -i $MATSIZE -r $MATSIZE -k $KVAL -o 1 >> ./results/n${MATSIZE}/results_1.txt
./build/krylov.exe -m import -n $MATSIZE -i $MATSIZE -r $MATSIZE -k $KVAL -o 1 >> ./results/n${MATSIZE}/results_1.txt
echo "" >> ./results/n${MATSIZE}/results_1.txt
done

MATSIZE=1024
for KVAL in {1..10}
do
echo ./build/krylov.exe -m import -n $MATSIZE -i $MATSIZE -r $MATSIZE -k $KVAL -o 2 >> ./results/n${MATSIZE}/results_2.txt
./build/krylov.exe -m import -n $MATSIZE -i $MATSIZE -r $MATSIZE -k $KVAL -o 2 >> ./results/n${MATSIZE}/results_2.txt
echo "" >> ./results/n${MATSIZE}/results_2.txt
done

MATSIZE=4096
for KVAL in {1..10}
do
echo ./build/krylov.exe -m import -n $MATSIZE -i $MATSIZE -r $MATSIZE -k $KVAL -o 2 >> ./results/n${MATSIZE}/results_2.txt
./build/krylov.exe -m import -n $MATSIZE -i $MATSIZE -r $MATSIZE -k $KVAL -o 2 >> ./results/n${MATSIZE}/results_2.txt
echo "" >> ./results/n${MATSIZE}/results_2.txt
done

MATSIZE=16384
for KVAL in {1..10}
do
echo ./build/krylov.exe -m import -n $MATSIZE -i $MATSIZE -r $MATSIZE -k $KVAL -o 2 >> ./results/n${MATSIZE}/results_2.txt
./build/krylov.exe -m import -n $MATSIZE -i $MATSIZE -r $MATSIZE -k $KVAL -o 2 >> ./results/n${MATSIZE}/results_2.txt
echo "" >> ./results/n${MATSIZE}/results_2.txt
done