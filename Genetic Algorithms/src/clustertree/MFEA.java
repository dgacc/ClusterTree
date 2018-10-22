package clustertree;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

import javax.swing.JFrame;

import dislay.Paint;
import filesinout.ReadFiles;
import operator.ChromosomeCmp;
import operator.Crossover;
import operator.Evaluation;
import operator.InitializeChromosome;
import operator.Mutations;
import random.MyRandom;
import structures.Cluster;
import structures.Individual;

public class MFEA {
	ReadFiles re = new ReadFiles();
	InitializeChromosome initChromo = new InitializeChromosome();
	Evaluation eva = new Evaluation();
	Individual chromo = new Individual();
	Crossover crossover = new Crossover();
	Mutations mutation = new Mutations();
	private static double crossOverRate = 0.8;
	public int maxGroup;
	public static int defaultPopLength = 100;
	ArrayList<Individual> population = new ArrayList<Individual>();

	public static void main(String[] args) {
		JFrame gf = new JFrame();
		gf.setVisible(true);
		gf.setSize(800, 800);
		gf.setTitle("the best Individual");
		gf.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		gf.setVisible(true);
//		ReadFiles.clusterReadFiles("C:/Users/TrungTB/Desktop/test/" +args[0]+ ".clt");
		ReadFiles.clusterReadFiles("C:/Users/TrungTB/Desktop/test/10eil76.clt");

		  PrintWriter pw = null;
			try {
				pw = new PrintWriter(new FileWriter(new File(args[1]), true));
			} catch (IOException e) {
				e.printStackTrace();
			}


		for(int loop = 0; loop < 1; loop ++){
		MyRandom.setSeed(loop);
		MFEA mainClass = new MFEA();
		int num_vertex = ReadFiles.num_vertex;
		int numberOfCluster = ReadFiles.numberOfCluster;
		ArrayList<Cluster> clusters = ReadFiles.clusters;
		ArrayList<Individual> pop = init();



		System.out.print(pop.get(0).getFactorialCost()[1] + " ");
		sortByCostIndex(pop, 0);
		System.out.println(pop.get(0).getFactorialCost()[0]);
		for (int i = 0; i < 7000; i++) {
			ArrayList<Individual> tempPop = new ArrayList<Individual>();
			while (tempPop.size() < defaultPopLength) {
				double r = 0 + (1 - 0) * MyRandom.r.nextDouble();
				int par1 = MyRandom.r.nextInt(defaultPopLength);
				int par2 = MyRandom.r.nextInt(defaultPopLength);
				if (pop.get(par1).getSkillFactor() == pop.get(par2).getSkillFactor() || r < crossOverRate) {
					Individual child = new Individual();
					child.setGene(mainClass.crossover.ClusterBFSCrossover(pop.get(par1).getGene(),
							pop.get(par1).getGene(), num_vertex,  clusters, MyRandom.r));
					tempPop.add(child);

				} else {
//					pop.get(par1)
//							.setGene(mainClass.mutation.mutationClusterTree(pop.get(par1).getGene(),
//									num_vertex, clusters,MyRandom.r ));
//					tempPop.add(pop.get(par1));
//					pop.get(par2)
//							.setGene(mainClass.mutation.mutationClusterTree(pop.get(par2).getGene(),
//									num_vertex, clusters,MyRandom.r));
//					tempPop.add(pop.get(par2));
				}
			}
			if (tempPop.size() > defaultPopLength) {
				tempPop.remove(MyRandom.r.nextInt(defaultPopLength));
			}
			calculate(tempPop);
			evaluatePopulation(tempPop);
			sortByCostIndex(tempPop, 0);

			Evaluation evaluation = new Evaluation();
			Paint  p = new Paint();
//			p.weightMatrix = tempPop.get(0).getGene();
//			p.fitness = evaluation.evaluation(ReadFiles.weightMatrix, tempPop.get(0).getGene(), ReadFiles.num_vertex, ReadFiles.root);
			p.num_vertex = num_vertex;
		    gf.add(p);
		    gf.setVisible(true);
			ArrayList<Individual> newPop = new ArrayList<Individual>();
			for (int j = 0; j < defaultPopLength / 2; j++) {
				newPop.add(pop.get(j));
				newPop.add(tempPop.get(j));
			}
			pop = newPop;



		}
		calculate(pop);
		evaluatePopulation(pop);

//		System.out.print(pop.get(0).getFactorialCost()[1] + " ");
		sortByCostIndex(pop, 0);
//		System.out.println(pop.get(0).getFactorialCost()[0]);
		Evaluation evaluation = new Evaluation();
//		 pw.println( evaluation.evaluation(ReadFiles.weightMatrix, pop.get(0).getGene(), ReadFiles.num_vertex, ReadFiles.root));
		  }
//		  pw.close();

	}

	public static void evaluatePopulation(ArrayList<Individual> pop) {
		sortByCostIndex(pop, 0);
		for (int i = 0; i < defaultPopLength; i++) {
			pop.get(i).setFactorialRank(i + 1, 0);
//			 System.out.println(pop.get(i).getFactorialCost()[0] + " ");
		}
		sortByCostIndex(pop, 1);
		for (int i = 0; i < defaultPopLength; i++) {
			pop.get(i).setFactorialRank(i + 1, 1);
		}
		for (int i = 0; i < defaultPopLength; i++) {
			if (pop.get(i).getFactorialRank()[0] < pop.get(i).getFactorialRank()[1]) {
				pop.get(i).setSkillFactor(0);
				pop.get(i).setScalarFitness(1.0 / (pop.get(i).getFactorialRank()[0]));
			} else {
				pop.get(i).setSkillFactor(1);
				pop.get(i).setScalarFitness(1.0 / (pop.get(i).getFactorialRank()[1]));
			}
		}
		Collections.sort(pop, ChromosomeCmp.compareByScalarFitness);
	}


	public static ArrayList<Individual> init() {
		MFEA mc = new MFEA();
		int numOfCity = ReadFiles.num_vertex;
		int numOfCluster = ReadFiles.numberOfCluster;
		double weightMatrix[][] = ReadFiles.weightMatrix;
		ArrayList<Cluster> clusters = ReadFiles.clusters;
		int maxGroupValue = 0;
		for (int i = 0; i < numOfCluster; i++) {
			if (maxGroupValue < clusters.get(i).getCluster().size()) {
				maxGroupValue = clusters.get(i).getCluster().size();
			}
		}
		mc.maxGroup = maxGroupValue;
		for (int i = 0; i < defaultPopLength; i++) {
			Individual c = new Individual();
			c.setGene(
					mc.initChromo.clusterPrimRST(weightMatrix, clusters, numOfCity));

			mc.population.add(c);
		}

		calculate(mc.population);
		return mc.population;
	}

	public static void calculate(ArrayList<Individual> pop) {
		MFEA mc = new MFEA();
		int num_vertex = ReadFiles.num_vertex;
		int numOfCluster = ReadFiles.clusters.size();
		double weightMatrix[][] = ReadFiles.weightMatrix;
	    ArrayList<Cluster> clusters = ReadFiles.clusters;
		int startVertex = ReadFiles.root;
		int maxGroupValue = 0;
		for (int i = 0; i < numOfCluster ; i++) {
			if (maxGroupValue < clusters.get(i).getCluster().size()) {
				maxGroupValue = clusters.get(i).getCluster().size();
			}
		}
		for (int i = 0; i < defaultPopLength; i++) {
			double[] temp = new double[2];
			temp[0] = mc.eva.evaluation( weightMatrix,pop.get(i).getGene(), num_vertex, startVertex);
			temp[1] = mc.eva.evaluation(
					weightMatrix, mc.eva.decodingMFOVertexInSubGraph(pop.get(i).getGene(), maxGroupValue, maxGroupValue),
					maxGroupValue, MyRandom.r.nextInt(maxGroupValue));
			pop.get(i).setFactorialCost(temp);
		}

	}

	public static void sortByCostIndex(ArrayList<Individual> pop, int index) {
		for (int i = 0; i < defaultPopLength; i++) {
			pop.get(i).cost = pop.get(i).getFactorialCost()[index];
		}
		Collections.sort(pop, ChromosomeCmp.compareByFactorialCost);

	}

	public static void printArray(double[][] arr) {
		for (int i = 0; i < arr.length; i++) {
			for (int j = 0; j < arr[i].length; j++) {
				System.out.print(arr[i][j] + " ");
			}
			System.out.println();
		}

	}
}
