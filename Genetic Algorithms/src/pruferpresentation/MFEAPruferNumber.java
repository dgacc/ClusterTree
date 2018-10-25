package pruferpresentation;

import java.io.File;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;
import java.util.concurrent.TimeUnit;

import dislay.Paint;
import dislay.Windows;
import edgepresentation.MEF2Instances;
import filesinout.ReadFiles;
import operator.ChromosomeCmp;
import operator.Crossover;
import operator.Evaluation;
import operator.InitializeChromosome;
import operator.Mutations;
import random.MyRandom;

//"10i30-17","10i30-17","10i45-18","10i45-18","10i60-21","10i60-21","10i65-21","10i65-21","10i70-21","10i70-21","10i75-22","10i75-22","10i90-33","10i90-33"};
//"4i200a","4i200h","4i200x1","4i200x2","4i200z","4i400a","4i400h","4i400x1","4i400x2","4i400z"
//"10berlin52","10eil51","10eil76","10kroB100","10pr76","10rat99","10st70","15berlin52" type 1
//"10berlin52","10eil51","10eil76","10kroB100","10pr76","10rat99",
//"15eil51","15eil76","15pr76","15st70","25eil101","25kroA100","25lin105","25rat99","50eil101","50kroA100",
//"50kroB100","50lin105","50rat99","5berlin52","5eil51","5eil76","5pr76","5st70"
import structures.Cluster;
import structures.Individual;

public class MFEAPruferNumber {
	ReadFiles re = new ReadFiles();
	private Evaluation eva = new Evaluation();
	private Mutations mutation = new Mutations();
	private Crossover crossover = new Crossover();
	private static InitializeChromosome initChromo = new InitializeChromosome();
	private static double rmp = 0.5;
	// private double muRate = 0.05;
	public int maxGroup;
	public static int defaultPopLength = 100;

	public static void main(String[] args) {
//		String str = args[0];
//		String[] test = str.split("_");
		 String[] test = { "15eil51", "15eil76" };
		 ReadFiles.clusterReadFiles("C:/Users/TrungTB/Desktop/test/" + test[0]
		 + ".clt");
		 ReadFiles.clusterReadFiles1("C:/Users/TrungTB/Desktop/test/" +
		 test[1] + ".clt");
//		ReadFiles.clusterReadFiles("test/" + test[0] + ".clt");
//		ReadFiles.clusterReadFiles1("test/" + test[1] + ".clt");
		MFEAPruferNumber mfeaPruferNumber = new MFEAPruferNumber();

		System.out.println(
				"------------------------------------------------------------------------------------------------------");
		// write file to use later
		String dirname = "Results/" + test[0] + "_and_" + test[1];
		File d = new File(dirname);
		// Bay gio tao thu muc.
		d.mkdirs();
		for (int loop = 0; loop < 30; loop++) {
			PrintWriter pw = null;
			PrintWriter pw1 = null;
			PrintWriter pw3 = null;
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

			// start MFEA algorithms here.
			ArrayList<Individual> currentPop = new ArrayList<Individual>();
			// Pre-processing cluster(sort ascending the number of elements)
			Collections.sort(ReadFiles.clusters, ChromosomeCmp.compareByNumberOfCluster);
			Collections.sort(ReadFiles.clusters1, ChromosomeCmp.compareByNumberOfCluster);
			// Initialize Information to new Individual in search space.l,
			initChromo.clutersVerticesInf(ReadFiles.clusters, ReadFiles.clusters1);
			initChromo.buildCluster();
			double muRate = mfeaPruferNumber.calMuRate(InitializeChromosome.maxClusters);
			// System.out.println(muRate);

			// initialize new population
			currentPop = mfeaPruferNumber.initiaizePopulationTest(InitializeChromosome.maxClusters, defaultPopLength,
					MyRandom.r);
			// Evaluate every Individual respects to every optimization task
			// in the multitasking environment.
			calculate(currentPop, defaultPopLength);
			bestIndividual = calculateFactorialRank(currentPop, defaultPopLength);
			calculateScalarFitness(currentPop, defaultPopLength);
			bestCost[0] = bestIndividual.get(0).getFactorialCost()[0];
			bestCost[1] = bestIndividual.get(1).getFactorialCost()[1];
			//
//			Windows windows = new Windows();
//			windows.runWindow(" cac cay");
			for (int i = 0; i < 500; i++) {
				ArrayList<Individual> offspringPop = new ArrayList<Individual>();
				ArrayList<Individual> intermiatePop = new ArrayList<Individual>();
				while (offspringPop.size() < defaultPopLength) {
					// Apply genetic operators on current-pop to generate an
					// offspring-pop (C).
					int par1 = MyRandom.r.nextInt(defaultPopLength);
					int par2 = MyRandom.r.nextInt(defaultPopLength);
					while (par1 == par2) {
						par2 = MyRandom.r.nextInt(defaultPopLength);
					}
					int geneLength = mfeaPruferNumber.calCulateGenLength(InitializeChromosome.maxClusters);
					offspringPop.addAll(mfeaPruferNumber.applyGeneticOperators(currentPop.get(par1),
							currentPop.get(par2), geneLength, muRate, MyRandom.r));

				}
				// Evaluate the individual in offspringPop for only selected
				// optimization tasks
				calculate(offspringPop, defaultPopLength);
				// . Concatenate offspring-pop and current-pop to form an
				// intermediate-pop (P ∪ C).
				intermiatePop.addAll(offspringPop);
				intermiatePop.addAll(currentPop);
				// iv. Update the scalar fitness (φ) and skill factor (τ) of
				// every individual in intermediate-pop.
				bestIndividual = (calculateFactorialRank(intermiatePop, defaultPopLength * 2));
				calculateScalarFitness(intermiatePop, defaultPopLength * 2);
				if (bestIndividual.get(0).getFactorialCost()[0] < bestCost[0]) {
					bestCost[0] = bestIndividual.get(0).getFactorialCost()[0];

				}
				pw3.print(i + "\t" + bestCost[0] + "\t");
				if (bestIndividual.get(1).getFactorialCost()[1] < bestCost[1]) {
					bestCost[1] = bestIndividual.get(1).getFactorialCost()[1];
				}
				// System.out.println(i +" \t" +bestCost[0]+ "\t "+bestCost[1]);
				pw3.println(bestCost[1]);

//				Paint p = new Paint();
//				double[][] caytest = mfeaPruferNumber.eva.decodingCayleyTest(bestIndividual.get(0).getGene1(),
//						mfeaPruferNumber.eva.getStartPoint(InitializeChromosome.maxClusters), ReadFiles.clusters,
//						InitializeChromosome.maxClusters, ReadFiles.num_vertex, MyRandom.r);
//				//
//				p.setPaint(caytest, ReadFiles.vertices, ReadFiles.clusters, ReadFiles.num_vertex,
//						bestIndividual.get(0).getFactorialCost()[0], 0, ReadFiles.root);
//				windows.addPaint(p);

				// v. Select the fittest individuals from intermediate-pop
				// to form the next current-pop (P)
				Collections.sort(intermiatePop, ChromosomeCmp.compareByScalarFitness);
				for (int j = 0; j < defaultPopLength; j++) {
					currentPop.set(j, intermiatePop.get(j));
				}

			}
			pw3.close();

			long end = System.currentTimeMillis();
			System.out.println("|" + loop + "\t|" + test[0] + "\t|" + test[1] + "\t|" + bestCost[0] + "\t|"
					+ bestCost[1] + "\t|" + getDateFromMillis(end - start) + "|");
			pw.println("Name: " + test[0]);
			pw.println("Seed: " + loop);
			pw.println("Fitness: " + bestCost[0]);
			pw.println("Time: " + getDateFromMillis(end - start));
			// for(int m = 0; m < ReadFiles.num_vertex; m ++){
			// pw.println();
			// for( int n = 0; n < ReadFiles.num_vertex; n++){
			// if( bestIndividual.get(0).getGene1()[m][n] > 0){
			// bestIndividual.get(0).getGene1()[m][n] = 1.0f; }
			// pw.print((int)bestIndividual.get(0).getGene1()[m][n]+"\t");
			//
			// }
			// }
			pw.close();
			pw1.println("Name: " + test[1]);
			pw1.println("Seed: " + loop);
			pw1.println("Fitness: " + bestCost[1]);
			pw1.println("Time: " + getDateFromMillis(end - start));
			// for(int m = 0; m < ReadFiles.num_vertex1; m ++){
			// pw1.println();
			// for( int n = 0; n < ReadFiles.num_vertex1; n++){
			// if( bestIndividual.get(0).getGene1()[m][n] > 0){
			// bestIndividual.get(0).getGene1()[m][n] = 1.0f; }
			// pw1.print((int)bestIndividual.get(0).getGene1()[m][n]+"\t");
			//
			// }
			// }
			pw1.close();

		}

	}

	public static ArrayList<Individual> calculateFactorialRank(ArrayList<Individual> pop, int pop_size) {
		ArrayList<Individual> bestIndividual = new ArrayList<>();
		sortByCostIndex(pop, pop_size, 0);
		bestIndividual.add(pop.get(0));
		for (int i = 0; i < pop_size; i++) {
			pop.get(i).setFactorialRank(i + 1, 0);
			// System.out.println(pop.get(i).getFactorialCost()[0] + " ");
		}
		sortByCostIndex(pop, pop_size, 1);
		bestIndividual.add(pop.get(0));
		for (int i = 0; i < pop_size; i++) {
			pop.get(i).setFactorialRank(i + 1, 1);
		}
		return bestIndividual;
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

	}

	/**
	 * calculate factorial cost
	 * 

	 * @param pop
	 */
	public static void calculate(ArrayList<Individual> pop, int pop_size) {
		MEF2Instances mc = new MEF2Instances();

		double weightMatrix[][] = ReadFiles.weightMatrix;
		double[][] weightMatrix1 = ReadFiles.weightMatrix1;
		int startVertex1 = ReadFiles.root;
		int startVertex2 = ReadFiles.root1;
		int num_vertex1 = ReadFiles.num_vertex;
		int num_vertex2 = ReadFiles.num_vertex1;

		int[] startPoint = mc.eva.getStartPoint(InitializeChromosome.maxClusters);

		for (int i = 0; i < pop_size; i++) {
			double[] temp = new double[2];
			double[][] TreeForEachInst = null;
			double[][] TreeForEachInst1 = null;
			if (pop.get(i).getSkillFactor() == 0) {
				TreeForEachInst = mc.eva.decodingCayleyTest(pop.get(i).getGene1(), startPoint, ReadFiles.clusters,
						InitializeChromosome.maxClusters, num_vertex1, MyRandom.r);
				temp[0] = mc.eva.evaluation(weightMatrix, TreeForEachInst, num_vertex1, startVertex1);
				temp[1] = Double.MAX_VALUE;
			} else if (pop.get(i).getSkillFactor() == 1) {
				TreeForEachInst1 = mc.eva.decodingCayleyTest(pop.get(i).getGene1(), startPoint, ReadFiles.clusters1,
						InitializeChromosome.maxClusters, num_vertex2, MyRandom.r);
				temp[1] = mc.eva.evaluation(weightMatrix1, TreeForEachInst1, num_vertex2, startVertex2);
				temp[0] = Double.MAX_VALUE;

			} else {
				TreeForEachInst = mc.eva.decodingCayleyTest(pop.get(i).getGene1(), startPoint, ReadFiles.clusters,
						InitializeChromosome.maxClusters, num_vertex1, MyRandom.r);

				temp[0] = mc.eva.evaluation(weightMatrix, TreeForEachInst, num_vertex1, startVertex1);
				TreeForEachInst1 = mc.eva.decodingCayleyTest(pop.get(i).getGene1(), startPoint, ReadFiles.clusters1,
						InitializeChromosome.maxClusters, num_vertex2, MyRandom.r);
				temp[1] = mc.eva.evaluation(weightMatrix1, TreeForEachInst1, num_vertex2, startVertex2);

			}


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

	public ArrayList<Individual> initiaizePopulationTest(ArrayList<Cluster> clusters, int popLength, Random rnd) {
		ArrayList<Individual> population = new ArrayList<Individual>();
		for (int m = 0; m < popLength; m++) {
			ArrayList<Integer> temp = new ArrayList<Integer>();
			ArrayList<Integer> temp1 = new ArrayList<Integer>();
			ArrayList<Integer> temp3 = new ArrayList<Integer>();
			temp3.removeAll(temp3);

			Individual individual = new Individual();

			int numberOfClusters = clusters.size();
			int[] presentVertex = new int[numberOfClusters];

			for (int i = 0; i < numberOfClusters; i++) {
				presentVertex[i] = rnd.nextInt(InitializeChromosome.maxClusterVertices[i]);
				temp.add(i, i);
				int numberClusterVertex = clusters.get(i).getCluster().size();
				if (numberClusterVertex == 1) {
					temp1.add(0, -1);// - 1 the number of cluster = 1;
				} else if (numberClusterVertex == 2) {
					temp1.add(0, -2); // -2 the number of cluster = 2;
				} else {
					for (int j = 0; j < numberClusterVertex; j++) {
						temp1.add(j, j);

					}
					Collections.shuffle(temp1, rnd);
					temp1 = deletePruferNumber(temp1, rnd);
				}
				temp3.addAll(temp1);
				temp1.removeAll(temp1);

			}
			Collections.shuffle(temp, rnd);
			temp = deletePruferNumber(temp, rnd);
			temp3.addAll(temp);
			temp.removeAll(temp);
			int num = temp3.size() + numberOfClusters;
			// System.out.println(num);
			int k = 0;
			for (int i = 0; i < num; i++) {
				if (i < temp3.size()) {
					individual.getGene1()[i] = temp3.get(i);
				} else {
					individual.getGene1()[i] = presentVertex[k];
					k++;
				}
			}
			individual.setSkillFactor(rnd.nextInt(2));
			population.add(individual);
		}
		return population;
	}

	/**
	 * initialize population
	 * 
	 * @param clusters
	 *            : the set of vertices
	 * @param popLength

	 * @param rnd
	 * @return
	 */


	public ArrayList<Individual> initiaizePopulation(ArrayList<Cluster> clusters, int popLength, Random rnd) {
		ArrayList<Individual> population = new ArrayList<Individual>();
		for (int m = 0; m < popLength; m++) {
			ArrayList<Integer> temp = new ArrayList<Integer>();
			ArrayList<Integer> temp1 = new ArrayList<Integer>();
			ArrayList<Integer> temp3 = new ArrayList<Integer>();
			temp3.removeAll(temp3);

			Individual individual = new Individual();

			int numberOfClusters = clusters.size();
			int[] presentVertex = new int[numberOfClusters];

			for (int i = 0; i < numberOfClusters; i++) {

				temp.add(i, i);
				int numberClusterVertex = clusters.get(i).getCluster().size();
				presentVertex[i] = rnd.nextInt(InitializeChromosome.maxClusterVertices[i]);
				for (int j = 0; j < numberClusterVertex; j++) {
					temp1.add(j, j);
				}
				Collections.shuffle(temp1, rnd);
				temp1 = deletePruferNumber(temp1, rnd);
				temp3.addAll(temp1);
				temp1.removeAll(temp1);

			}
			Collections.shuffle(temp, rnd);
			temp = deletePruferNumber(temp, rnd);
			temp3.addAll(temp);
			temp.removeAll(temp);
			int num = temp3.size() + numberOfClusters;

			int k = 0;
			for (int i = 0; i < num; i++) {
				if (i < temp3.size()) {
					individual.getGene1()[i] = temp3.get(i);
				} else {
					individual.getGene1()[i] = presentVertex[k];
					k++;
				}
			}
			individual.setSkillFactor(rnd.nextInt(2));
			population.add(individual);
		}
		return population;
	}

	/**
	 * delete randomly from the number at first code. delete two last elements
	 * 
	 * @param individual
	 *            : the number at first
	 * @param rnd
	 *            : random seed
	 * @return the prufer number code after delete
	 */

	public ArrayList<Integer> deletePruferNumber(ArrayList<Integer> individual, Random rnd) {

		int geneLength = individual.size();
		individual.remove(geneLength - 1);
		individual.remove(geneLength - 2);
		return individual;
	}


	/**
	 * calculate the length of prufernumber gene;
	 * 
	 * @param maxClusters
	 * @return
	 */
	public int calCulateGenLength(ArrayList<Cluster> maxClusters) {
		int numberOfCluster = maxClusters.size();
		int geneLength = 0;
		for (int i = 0; i < numberOfCluster; i++) {
			int numberClusterVertices = maxClusters.get(i).getCluster().size();
			if (numberClusterVertices < 3) {
				geneLength++;
			} else {
				geneLength += (numberClusterVertices - 2);
			}
		}

		geneLength += (2 * numberOfCluster - 2);

		return geneLength;
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

	public ArrayList<Individual> applyGeneticOperators(Individual par1, Individual par2, int GeneLength, double muRate,
			Random rnd) {

		double r = rnd.nextDouble();
		ArrayList<Individual> Offsprings = new ArrayList<>();

		if (par1.getSkillFactor() == par2.getSkillFactor() || r < rmp) {

			Individual child = new Individual();
			Individual child1 = new Individual();
			child.setGene(
					crossover.pruferNumberCrossover(par1.getGene1(), par2.getGene1(), GeneLength, MyRandom.r).get(0));
			child1.setGene(
					crossover.pruferNumberCrossover(par1.getGene1(), par2.getGene1(), GeneLength, MyRandom.r).get(1));
			if (rnd.nextDouble() < 0.5) {
				child.setSkillFactor(par1.getSkillFactor());
			} else {
				child.setSkillFactor(par2.getSkillFactor());
			}
			if (rnd.nextDouble() < 0.5) {
				child1.setSkillFactor(par1.getSkillFactor());
			} else {
				child1.setSkillFactor(par2.getSkillFactor());
			}

			Offsprings.add(child);
			Offsprings.add(par2);

		} else {
			par1.setGene(mutation.mutationPrufer(InitializeChromosome.maxClusters, rnd, par1.getGene1(),
					eva.getStartPoint(InitializeChromosome.maxClusters), muRate));
			par2.setGene(mutation.mutationPrufer(InitializeChromosome.maxClusters, rnd, par2.getGene1(),
					eva.getStartPoint(InitializeChromosome.maxClusters), muRate));
			Offsprings.add(par1);
			Offsprings.add(par2);
		}
		return Offsprings;
	}

}
