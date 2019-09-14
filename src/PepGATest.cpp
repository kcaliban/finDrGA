/* Copyright 2019 Fabian Krause */
#include <gtest/gtest.h>

/**** Tests for GenAlgInst ****/
#include <random>
#include "lib/GenAlgInst.h"
#include "lib/Genome.h"
/** TEST #1: Top individuals being copied **/
class TestGenomeCpy : public Genome<int> {
 public:
    int crossOver(int& x, int& y, ...) {
      return -1;
    }
    int mutate(int& x, ...) {
      return -1;
    }
};

class TestFitnessFunctionCpy : public FitnessFunction<int> {
 public:
   float calculateFitness(int x, ...) {
     return x;
   }
};

TEST(GenAlgInst, topIndividuals50) {
  std::random_device rd;
  std::mt19937 mt(rd());
  TestFitnessFunctionCpy testFitnessFunction;
  TestGenomeCpy testGenome;
  GenAlgInst<int, TestGenomeCpy, TestFitnessFunctionCpy> genAlgInst(&mt);

  std::vector<int> initialPop = {10, 10, 5, 3};
  std::vector<int> nextgen = genAlgInst.nextGen(testGenome,
                                    testFitnessFunction,
                                    initialPop,
                                    0,
                                    0.5);
  ASSERT_EQ(2, std::count(nextgen.begin(), nextgen.end(), 10));
}

TEST(GenAlgInst, topIndividuals25) {
  std::random_device rd;
  std::mt19937 mt(rd());
  TestFitnessFunctionCpy testFitnessFunction;
  TestGenomeCpy testGenome;
  GenAlgInst<int, TestGenomeCpy, TestFitnessFunctionCpy> genAlgInst(&mt);

  std::vector<int> initialPop = {33, 15, 20, 39, 40, 50, 45, 69};
  std::vector<int> nextgen = genAlgInst.nextGen(testGenome,
                                    testFitnessFunction,
                                    initialPop,
                                    0,
                                    0.25);
  ASSERT_EQ(1, std::count(nextgen.begin(), nextgen.end(), 69));
  ASSERT_EQ(1, std::count(nextgen.begin(), nextgen.end(), 50));
}

TEST(GenAlgInst, topIndividuals0) {
  std::random_device rd;
  std::mt19937 mt(rd());
  TestFitnessFunctionCpy testFitnessFunction;
  TestGenomeCpy testGenome;
  GenAlgInst<int, TestGenomeCpy, TestFitnessFunctionCpy> genAlgInst(&mt);

  std::vector<int> initialPop = {33, 15, 20, 39, 40, 50, 45, 69};
  std::vector<int> nextgen = genAlgInst.nextGen(testGenome,
                                    testFitnessFunction,
                                    initialPop,
                                    0,
                                    0);
  ASSERT_EQ(8, std::count(nextgen.begin(), nextgen.end(), -1));
}
/** TEST #2: Probability of being recombined is roughly correct **/
class TestGenomeRec : public Genome<int> {
 public:
    int crossOver(int& x, int& y, ...) {
      return (x < y) ? y : x;
    }
    int mutate(int& x, ...) {
      return 0;
    }
};

class TestFitnessFunctionRec : public FitnessFunction<int> {
 public:
   float calculateFitness(int x, ...) {
     return x;
   }
};

TEST(GenAlgInst, recomb1) {
  std::random_device rd;
  std::mt19937 mt(rd());
  TestFitnessFunctionRec testFitnessFunction;
  TestGenomeRec testGenome;
  GenAlgInst<int, TestGenomeRec, TestFitnessFunctionRec> genAlgInst(&mt);

  std::vector<int> initialPop = {1, 1, 1, 10};
  int numberAllTens = 0;
  int n = 500;
  for (int i = 0; i < n; i++) {
    std::vector<int> nextgen = genAlgInst.nextGen(testGenome,
                                      testFitnessFunction,
                                      initialPop,
                                      0,
                                      0);
    if (std::count(nextgen.begin(), nextgen.end(), 10) == 4) {
      numberAllTens += 1;
    }
  }
  // Recombination: Taking the max
  // For one recombined individual we have:
  //   P(¬10)    = (3/13)^2
  //   P(10)     = (10/13) * (3/13) + (3/13) * (10/13) + (10/13)^2
  //             = 1 - (3/13)^2
  // With following probability we have all 10s in next generation:
  //   P(ALL_10) = 1 - sum_{i = 1}^4 (4 choose i) * P(¬10)^2i * P(10)^(4 - i)
  //             ≈ 0.8
  // Ergo roughly 80% of all generations should have all tens
  ASSERT_NEAR(0.8, (float) numberAllTens / (float) n, 0.05);
}

/** TEST #3: Probability of mutation is correct **/
class TestGenomeMut : public Genome<int> {
 public:
    int crossOver(int& x, int& y, ...) {
      return x;
    }
    int mutate(int& x, ...) {
      return -1;
    }
};

class TestFitnessFunctionMut : public FitnessFunction<int> {
 public:
   float calculateFitness(int x, ...) {
     return x;
   }
};

TEST(GenAlgInst, mutate50) {
  std::random_device rd;
  std::mt19937 mt(rd());
  TestFitnessFunctionMut testFitnessFunction;
  TestGenomeMut testGenome;
  GenAlgInst<int, TestGenomeMut, TestFitnessFunctionMut> genAlgInst(&mt);

  std::vector<int> initialPop = {10, 10, 10, 10, 10, 10, 10, 10};
  int numberMut = 0;
  int n = 500;
  for (int i = 0; i < n; i++) {
    std::vector<int> nextgen = genAlgInst.nextGen(testGenome,
                                      testFitnessFunction,
                                      initialPop,
                                      0.5,
                                      0);
    // Each individual has probability 1/2 to mutate
    // Expected number of mutations per generation is therefore 8 * 1/2 = 4
    numberMut += std::count(nextgen.begin(), nextgen.end(), -1);
  }
  // Number of mutations per generation averaged over number of generations
  // should be close to 4
  ASSERT_NEAR(4, (float) numberMut / (float) n, 0.1);
}

TEST(GenAlgInst, mutate25) {
  std::random_device rd;
  std::mt19937 mt(rd());
  TestFitnessFunctionMut testFitnessFunction;
  TestGenomeMut testGenome;
  GenAlgInst<int, TestGenomeMut, TestFitnessFunctionMut> genAlgInst(&mt);

  std::vector<int> initialPop = {10, 10, 10, 10, 10, 10, 10, 10};
  int numberMut = 0;
  int n = 500;
  for (int i = 0; i < n; i++) {
    std::vector<int> nextgen = genAlgInst.nextGen(testGenome,
                                      testFitnessFunction,
                                      initialPop,
                                      0.25,
                                      0);
    // Each individual has probability 1/2 to mutate
    // Expected number of mutations per generation is therefore 8 * 1/4 = 2
    numberMut += std::count(nextgen.begin(), nextgen.end(), -1);
  }
  // Number of mutations per generation averaged over number of generations
  // should be close to 2
  ASSERT_NEAR(2, (float) numberMut / (float) n, 0.1);
}

/**** PDB to FASTA tests ****/
#include "PoolManager/PoolManager.h"
#include "Info.h"

TEST(PoolMGR, PDBtoFASTA) {
  std::vector<std::string> blub;
  std::vector<std::pair<std::string, std::string>> PDBs;
  Info * info = new Info(false, false, "");
  PoolMGR poolmgr("", "", "", "", "", blub, 0, 0, "", "", "", "", "", "", 0, 0, info);
  PDBs.push_back(std::make_pair("SEAHTLLYGT", "src/testpdbs/10_5zop_D_6.pdb"));
  PDBs.push_back(std::make_pair("SLYPAQMYSR", "src/testpdbs/10_pdb2pyj_D_25.pdb"));
  PDBs.push_back(std::make_pair("PKFFRYISERY", "src/testpdbs/11_2y35_D_1.pdb"));
  PDBs.push_back(std::make_pair("SFYKNFTIDCS", "src/testpdbs/11_pdb1pd0_D_15.pdb"));
  PDBs.push_back(std::make_pair("GVLGKELIDVLD", "src/testpdbs/12_2r8j_D_17.pdb"));
  PDBs.push_back(std::make_pair("MHESMKLFDSIC", "src/testpdbs/12_pdb2g83_D_14.pdb"));
  PDBs.push_back(std::make_pair("GDGVEEAF", "src/testpdbs/8_5icn_D_18.pdb"));
  PDBs.push_back(std::make_pair("LVLLLRVL", "src/testpdbs/8_5vo8_D_14.pdb"));
  PDBs.push_back(std::make_pair("LPDLFEEQ", "src/testpdbs/8_5yko_D_16.pdb"));
  PDBs.push_back(std::make_pair("WDEHDFYCR", "src/testpdbs/9_5usn_D_2.pdb"));
  // PDBs.push_back(std::make_pair("
  for (std::vector<std::pair<std::string, std::string>>::const_iterator it = PDBs.begin();
      it != PDBs.end(); it++) {
    std::string fasta = it->first;
    std::string pdb = it->second;
    std::string out = poolmgr.PDBtoFASTA(it->second);
    ASSERT_EQ(fasta, out);
  }
  free(info);
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
