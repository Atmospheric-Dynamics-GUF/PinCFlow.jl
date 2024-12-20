import subprocess

#dir_1 = "/PF/runs/redo_tests/"
dir_1 = "/PF/runs/"
dir_2 = "/PF/runs_reftests/"
#dir_2 = "/PF/runs/runs_test2"

# Get host and user name.
host_name = subprocess.getoutput("hostname")
user_name = subprocess.getoutput("whoami")

if "levante" in host_name:
  # Levante cluster
  data_path = "/scratch/b/" + user_name + dir_1
  reference_path = "/scratch/b/" + user_name + dir_2

elif "login" in host_name:
  # Goethe cluster
  data_path = "/scratch/atmodynamics/" + user_name + dir_1
  reference_path = "/scratch/atmodynamics/" + user_name + dir_2

else:
  # Local machine
  data_path = "/home/" + user_name + dir_1
  reference_path = "/home/" + user_name + dir_2
  
