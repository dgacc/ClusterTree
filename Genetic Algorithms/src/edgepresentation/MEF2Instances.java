package edgepresentation;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.concurrent.TimeUnit;

import dislay.Windows;
import filesinout.ReadFiles;
import operator.ChromosomeCmp;
import operator.Crossover;
import operator.Evaluation;
import operator.InitializeChromosome;
import operator.Mutations;
import random.MyRandom;
import structures.Cluster;
import structures.Individual;

public class MEF2Instances {
	ReadFiles re = new ReadFiles();
	public Evaluation eva = new Evaluation();
	protected Individual chromo = new Individual();
	protected Mutations mutation = new Mutations();
	protected Crossover crossover = new Crossover();
	protected static double crossOverRate = 0.5;
	protected int maxGroup;
	protected static int defaultPopLength = 100;
	protected ArrayList<Individual> population = new ArrayList<Individual>();

	public static void main(String[] args) {

		// Tao window cho bai toan
		// Windows windows = new Windows();

		// windows.runWindow(" cac cay");
		// tot nhat trong cac the he
		String str = args[0];
		String[] test = str.split("-");
		// String[] test = {"2lin105-2x1"};
		// for(int k= 0; k < test.length; k= k+2){
		// ReadFiles.clusterReadFiles("C:/Users/TrungTB/Desktop/test/"+
		// test[k]+".clt");
		ReadFiles.clusterReadFiles("test/" + test[0] + ".clt");
		ReadFiles.clusterReadFiles1("test/" + test[1] + ".clt");

		System.out.println("-------------------------------------------------------------------------------");
		// write file to use later
		String dirname = "Results/" + test[0] + "_and_" + test[1];
		File d = new File(dirname);
		// Bay gio tao thu muc.
		d.mkdirs();
		for (int loop = 0; loop < 30; loop++) {
			PrintWriter pw = null;
			PrintWriter pw1 = null;
			PrintWriter pw3 = null;
			// PrintWriter pw4 = null;

			try {
				pw = new PrintWriter(new FileWriter(dirname + "/" + test[0] + "-seed(" + loop + ").opt", true));
				pw1 = new PrintWriter(new FileWriter(dirname + "/" + test[1] + "-seed(" + loop + ").opt", true));
				pw3 = new PrintWriter(
						new FileWriter(dirname + "/" + test[0] + "_and_" + test[1] + "_seed(" + loop + ").gen", true));
			} catch (IOException e) {
				e.printStackTrace();
			}
			pw3.println("Generations \t" + test[0] + "\t\t" + test[1]);
			long start = System.currentTimeMillis();
			double[] bestCost = new double[2];
			ArrayList<Individual> bestIndividual = new ArrayList<Individual>();
			MyRandom.setSeed(loop);
			MEF2Instances mef2 = new MEF2Instances();
			InitializeChromosome initChromo = new InitializeChromosome();
			initChromo.clutersVerticesInf(ReadFiles.clusters, ReadFiles.clusters1);
			initChromo.buildCluster();
			// double muRate =mef2.calMuRate(InitializeChromosome.maxClusters);
			double muRate = 0.05;
			ArrayList<Individual> pop = new ArrayList<Individual>();
			Individual individual = new Individual();
			// ininitlize new population
			for (int i = 0; i < defaultPopLength; i++) {
				individual.setGene(initChromo.initializeGroupGene(InitializeChromosome.maxClusters, MyRandom.r));

				pop.add(individual);

			}

			calculate(pop, defaultPopLength);
			calculateFactorialRank(pop, defaultPopLength);
			calculateScalarFitness(pop, defaultPopLength);

			sortByCostIndex(pop, defaultPopLength, 0);
			bestCost[0] = pop.get(0).getFactorialCost()[0];
			// System.out.print(pop.get(0).getFactorialCost()[0] + " ");
			bestIndividual.add(0, pop.get(0));
			sortByCostIndex(pop, defaultPopLength, 1);
			bestCost[1] = pop.get(0).getFactorialCost()[1];
			// System.out.println(pop.get(0).getFactorialCost()[1]);
			bestIndividual.add(1, pop.get(0));

			Windows windows1 = new Windows();

			for (int i = 0; i < 500; i++) {
				ArrayList<Individual> newPop = new ArrayList<Individual>();
				// Windows windows = new Windows();
				// windows.runWindow(" cac cay" + i);
				// System.out.println(i);
				ArrayList<Individual> tempPop = new ArrayList<Individual>();
				while (tempPop.size() < defaultPopLength) {
					double r = MyRandom.r.nextDouble();
					int par1 = MyRandom.r.nextInt(defaultPopLength);
					int par2 = MyRandom.r.nextInt(defaultPopLength);
					if (pop.get(par1).getSkillFactor() == pop.get(par2).getSkillFactor() || r < crossOverRate) {
						Individual child = new Individual();
						child.setGene(mef2.crossover.clusterCrossover1(pop.get(par1).getGene(), pop.get(par2).getGene(),
								InitializeChromosome.maxNumberVertices, InitializeChromosome.maxClusters,
								ReadFiles.clusters.size(), ReadFiles.clusters1.size(),
								InitializeChromosome.minClusterVertices, InitializeChromosome.maxClusterVertices,
								MyRandom.r, windows1));
						if (r < 0.5) {
							child.setSkillFactor(pop.get(par1).getSkillFactor());
						} else {
							child.setSkillFactor(pop.get(par2).getSkillFactor());
						}
						tempPop.add(child);

					} else {

						pop.get(par1)
								.setGene(mef2.mutation.mutationClusterTree(pop.get(par1).getGene(),
										InitializeChromosome.maxNumberVertices, InitializeChromosome.maxClusters,
										InitializeChromosome.minClusterVertices,
										InitializeChromosome.maxClusterVertices, MyRandom.r, muRate));
						tempPop.add(pop.get(par1));
						pop.get(par2)
								.setGene(mef2.mutation.mutationClusterTree(pop.get(par2).getGene(),
										InitializeChromosome.maxNumberVertices, InitializeChromosome.maxClusters,
										InitializeChromosome.minClusterVertices,
										InitializeChromosome.maxClusterVertices, MyRandom.r, muRate));
						tempPop.add(pop.get(par2));
					}
				}
				if (tempPop.size() > defaultPopLength) {
					tempPop.remove(defaultPopLength);
				}
				calculate(tempPop, defaultPopLength);

				newPop.addAll(tempPop);
				newPop.addAll(pop);

				calculateFactorialRank(newPop, defaultPopLength * 2);

				sortByCostIndex(newPop, defaultPopLength * 2, 0);
				if (newPop.get(0).getFactorialCost()[0] < bestCost[0]) {
					bestCost[0] = newPop.get(0).getFactorialCost()[0];
					bestIndividual.set(0, newPop.get(0));
				}
				pw3.print(i + "\t" + bestCost[0] + "\t");
				// System.out.print(newPop.get(0).getFactorialCost()[0] + " ");
				sortByCostIndex(newPop, defaultPopLength * 2, 1);
				if (newPop.get(0).getFactorialCost()[1] < bestCost[1]) {
					bestCost[1] = newPop.get(0).getFactorialCost()[1];
					bestIndividual.set(1, newPop.get(0));
				}
				pw3.println(bestCost[1]);

				// System.out.println(newPop.get(0).getFactorialCost()[1]);

				// Paint p = new Paint();
				// ArrayList<double[][]> caytest =
				// initChromo.decodingTwoTree(newPop.get(0).getGene(),
				// ReadFiles.clusters, ReadFiles.clusters1,
				// ReadFiles.num_vertex, ReadFiles.num_vertex1, MyRandom.r);
				//
				// p.setPaint(caytest.get(1), ReadFiles.vertices1,
				// ReadFiles.clusters1, ReadFiles.num_vertex1,
				// tempPop.get(0).getFactorialCost()[1], 0, ReadFiles.root1);
				// windows.addPaint(p);
				//

				calculateScalarFitness(newPop, defaultPopLength * 2);

				for (int j = 0; j < defaultPopLength; j++) {
					pop.set(j, newPop.get(j));

				}

			}
			pw3.close();

			long end = System.currentTimeMillis();
			// System.out.print(bestCost[0] + " ");
			// System.out.println();
			System.out.println("|" + loop + "\t|" + test[0] + "\t|" + test[1] + "\t|" + bestCost[0] + "\t|"
					+ bestCost[1] + "\t|" + getDateFromMillis(end - start) + "|");
			pw.println("Name: " + test[0]);
			pw.println("Seed: " + loop);
			pw.println("Fitness: " + bestCost[0]);
			pw.println("Time: " + getDateFromMillis(end - start));
			for (int m = 0; m < ReadFiles.num_vertex; m++) {
				pw.println();
				for (int n = 0; n < ReadFiles.num_vertex; n++) {
					if (bestIndividual.get(0).getGene()[m][n] > 0) {
						bestIndividual.get(0).getGene()[m][n] = 1.0f;
					}
					pw.print((int) bestIndividual.get(0).getGene()[m][n] + "\t");

				}
			}
			pw.close();

			pw1.println("Name: " + test[1]);
			pw1.println("Seed: " + loop);
			pw1.println("Fitness: " + bestCost[1]);
			pw1.println("Time: " + getDateFromMillis(end - start));
			for (int m = 0; m < ReadFiles.num_vertex1; m++) {
				pw1.println();
				for (int n = 0; n < ReadFiles.num_vertex1; n++) {
					if (bestIndividual.get(0).getGene()[m][n] > 0) {
						bestIndividual.get(0).getGene()[m][n] = 1.0f;
					}
					pw1.print((int) bestIndividual.get(0).getGene()[m][n] + "\t");

				}
			}
			pw1.close();

		}

	}
	// }

	public static void calculateFactorialRank(ArrayList<Individual> pop, int pop_size) {
		sortByCostIndex(pop, pop_size, 0);
		for (int i = 0; i < pop_size; i++) {
			pop.get(i).setFactorialRank(i + 1, 0);
			// System.out.println(pop.get(i).getFactorialCost()[0] + " ");
		}
		sortByCostIndex(pop, pop_size, 1);
		for (int i = 0; i < pop_size; i++) {
			pop.get(i).setFactorialRank(i + 1, 1);
		}
	}

	public static void calculateScalarFitness(ArrayList<Individual> pop, int pop_size) {
		for (int i = 0; i < pop_size; i++) {
			if (pop.get(i).getFactorialRank()[0] < pop.get(i).getFactorialRank()[1]) {
				pop.get(i).setSkillFactor(0);
				pop.get(i).setScalarFitness(1.0 / (pop.get(i).getFactorialRank()[0]) + 1);
			} else {
				pop.get(i).setSkillFactor(1);
				pop.get(i).setScalarFitness(1.0 / (pop.get(i).getFactorialRank()[1]) + 1);
			}
		}
		Collections.sort(pop, ChromosomeCmp.compareByScalarFitness);
	}

	/**
	 * calculate factorial cost
	 * 
	 * @param pop
	 */
	public static void calculate(ArrayList<Individual> pop, int pop_size) {
		MEF2Instances mc = new MEF2Instances();
		InitializeChromosome initChromo = new InitializeChromosome();
		double weightMatrix[][] = ReadFiles.weightMatrix;
		double[][] weightMatrix1 = ReadFiles.weightMatrix1;
		int startVertex1 = ReadFiles.root;
		int startVertex2 = ReadFiles.root1;
		int num_vertex1 = ReadFiles.num_vertex;
		int num_vertex2 = ReadFiles.num_vertex1;

		for (int i = 0; i < pop_size; i++) {
			double[] temp = new double[2];
			ArrayList<double[][]> TreeForEachInst = new ArrayList<double[][]>();
			TreeForEachInst = initChromo.decodingTwoTree(pop.get(i).getGene(), ReadFiles.clusters, ReadFiles.clusters1,
					num_vertex1, num_vertex2, MyRandom.r);

			temp[0] = mc.eva.evaluation(weightMatrix, TreeForEachInst.get(0), num_vertex1, startVertex1);
			temp[1] = mc.eva.evaluation(weightMatrix1, TreeForEachInst.get(1), num_vertex2, startVertex2);
			pop.get(i).setFactorialCost(temp);
		}

	}

	public static void sortByCostIndex(ArrayList<Individual> pop, int pop_size, int index) {
		for (int i = 0; i < pop_size; i++) {
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

	public static String getDateFromMillis(long millis) {
		String string = String.format("%02d:%02d:%02d.%03d", TimeUnit.MILLISECONDS.toHours(millis),
				TimeUnit.MILLISECONDS.toMinutes(millis)
						- TimeUnit.HOURS.toMinutes(TimeUnit.MILLISECONDS.toHours(millis)),
				TimeUnit.MILLISECONDS.toSeconds(millis)
						- TimeUnit.MINUTES.toSeconds(TimeUnit.MILLISECONDS.toMinutes(millis)),
				millis - TimeUnit.SECONDS.toMillis(TimeUnit.MILLISECONDS.toSeconds(millis)));
		return string;
	}

	public double calMuRate(ArrayList<Cluster> maxClusters) {

		int numberOfCluster = maxClusters.size();
		int numberVertices = 0;
		for (int i = 0; i < numberOfCluster; i++) {
			int numberClusterVertices = maxClusters.get(i).getCluster().size();
			numberVertices += (numberClusterVertices);
		}
		return 1.0 / numberVertices;
	}

}