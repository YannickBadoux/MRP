start_time=$(date +%s)

python /home/ybadoux/Documents/MRP/MRP_code/acs_star_planet_moon.py --a_sp 1 --density 3 --output /home/ybadoux/Documents/MRP/acs_spm_test
python /home/ybadoux/Documents/MRP/MRP_code/acs_star_planet_moon.py --a_sp 5 --density 3 --output /home/ybadoux/Documents/MRP/acs_spm_test
python /home/ybadoux/Documents/MRP/MRP_code/acs_star_planet_moon.py --a_sp 10 --density 3 --output /home/ybadoux/Documents/MRP/acs_spm_test
python /home/ybadoux/Documents/MRP/MRP_code/acs_star_planet_moon.py --a_sp 50 --density 3 --output /home/ybadoux/Documents/MRP/acs_spm_test
python /home/ybadoux/Documents/MRP/MRP_code/acs_star_planet_moon.py --a_sp 100 --density 3 --output /home/ybadoux/Documents/MRP/acs_spm_test

end_time=$(date +%s)
runtime=$((end_time - start_time))
echo "Total runtime: $runtime seconds"