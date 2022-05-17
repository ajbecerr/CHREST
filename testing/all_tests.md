### To run all tests with CCR Academic Cluster Shell Access:

1. Log into CCR OnDemand

2. Access the Academic Cluster Shell

3. Clone this Repository (if needed)
  ```bash
	git clone https://github.com/ajbecerr/CHREST.git
  ```
  
4. Execute all tests and save a copy of all ouput(s)
  ```bash
	cd CHREST/testing/
	bash all_tests.sh |& tee all_tests_Academic_Cluster_Shell_Access.out
  ```

### To run all tests with an Interactive Session:

1. Log into CCR OnDemand

2. Request an Interative Session (Partition/QOS = general-compute, Node Features = CPU-Gold-6130)

3. Clone this Repository (if needed)
  ```bash
	git clone https://github.com/ajbecerr/CHREST.git
  ```
  
4. Execute all tests and save a copy of all ouput(s)
  ```bash
	cd CHREST/testing/
	bash all_tests.sh |& tee all_tests_Academic_Cluster_Interactive_Session.out
  ```