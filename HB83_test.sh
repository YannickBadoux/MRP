start_time=$(date +%s)

python HB83.py --velocity 1 --n_b 100 --sim_per_b 100 --output /home/ybadoux/Documents/MRP/HB83_test2
python HB83.py --velocity 1.16 --n_b 100 --sim_per_b 100 --output /home/ybadoux/Documents/MRP/HB83_test2
python HB83.py --velocity 1.34 --n_b 100 --sim_per_b 100 --output /home/ybadoux/Documents/MRP/HB83_test2
python HB83.py --velocity 1.55 --n_b 100 --sim_per_b 100 --output /home/ybadoux/Documents/MRP/HB83_test2
python HB83.py --velocity 1.79 --n_b 100 --sim_per_b 100 --output /home/ybadoux/Documents/MRP/HB83_test2
python HB83.py --velocity 2.07 --n_b 100 --sim_per_b 100 --output /home/ybadoux/Documents/MRP/HB83_test2
python HB83.py --velocity 2.40 --n_b 100 --sim_per_b 100 --output /home/ybadoux/Documents/MRP/HB83_test2
python HB83.py --velocity 2.78 --n_b 100 --sim_per_b 100 --output /home/ybadoux/Documents/MRP/HB83_test2
python HB83.py --velocity 3.21 --n_b 100 --sim_per_b 100 --output /home/ybadoux/Documents/MRP/HB83_test2
python HB83.py --velocity 3.72 --n_b 100 --sim_per_b 100 --output /home/ybadoux/Documents/MRP/HB83_test2

end_time=$(date +%s)
runtime=$((end_time - start_time))
echo "Total runtime: $runtime seconds"