if [ -z "$1" ] || [ -z "$2" ]; then
    echo -e "usage ./run <test or ecoli> <compile?>\n"
    echo -e "\tEX: ./run test yes => compiles code and runs on test_data"
    exit 1
fi



if [$2 == "yes"]; then
    echo "compiling"
    g++ build_graphs.cpp -o build_graphs
    if [ $? -ne 0 ]; then
        echo "g++ encountered an error."
    exit 1
    fi
fi


if ["$1" == "test"]; then
    ./build_graphs ecoli_data/SRR10971019_sub.fasta ecoli_data/sl_raw_align.txt
else
    ./build_graphs test_data/pac.fasta test_data/sl_raw_align.txt
fi