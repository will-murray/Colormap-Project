if [ -z "$1" ] || [ -z "$2" ]; then
    echo -e "usage ./run <test or ecoli> <compile?>\n"
    echo -e "\tEX: ./run test yes => compiles code and runs on test_data"
    exit 1
fi



if [ "$2" == "yes" ]; then
    echo "Compiling"
    g++ colormap.cpp -o colormap
    if [ $? -ne 0 ]; then
        echo "g++ encountered an error."
        exit 1
    else
        echo -e "compilation success >:)"
    fi
fi


if [ "$1" == "test" ]; then
    echo "running ColoRMap on test_data/"
    ./colormap test_data/pac.fasta test_data/sl_raw_align.txt
else
    echo "running ColoRMap on ecoli_data/"
    ./colormap ecoli_data/SRR10971019_sub.fasta ecoli_data/sl_raw_align.txt
fi


echo "Updating graph.png"
dot -Tpng graph.dot -o graph.png