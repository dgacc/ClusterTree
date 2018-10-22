package pruferpresentation;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;

import edgepresentation.MEF2Instances;
import filesinout.ReadFiles;
import operator.ChromosomeCmp;
import operator.Crossover;
import operator.Evaluation;
import operator.InitializeChromosome;
import operator.Mutations;
import operator.Selection;
import random.MyRandom;
import structures.Individual;
import structures.Population;

public class GAPruferNumber {
	private static double m_rate = 0.05;
	private static final double p_rate = 0.5; // crossover rate
	private static int generation = 500;
	private static int populationLength = 100; // default length of population
	private static Crossover crossover = new Crossover();
	private static Mutations mutation = new Mutations();
	private static Selection selection = new Selection();
	private static Evaluation evaluation = new Evaluation();

	public static void main(String[] args) {
		MFEAPruferNumber mfeaPruferNumber = new MFEAPruferNumber();
		ReadFiles.clusterReadFiles("test/" + args[0] + ".clt");
//		 Windows windows = new Windows();
//		 windows.runWindow("The best ever");
		String dirname = "Result/" + args[0];
		File d1 = new File(dirname);
		// Bay gio tao thu muc.
		d1.mkdirs();

		//
		System.out.println("-------------------------------------------------------------------------------");

		for (int loop = 0; loop < 30; loop++) {
			MyRandom.setSeed(loop);
			// write to files
			PrintWriter pw = null;
			PrintWriter pw3 = null;

			try {
				pw = new PrintWriter(new FileWriter(dirname + "/" + args[0] + "-seed(" + loop + ").opt", true));
				pw3 = new PrintWriter(new FileWriter(dirname + "/" + args[0] + "_seed(" + loop + ").gen", true));
			} catch (IOException e) {
				e.printStackTrace();
			}
			pw3.println("Generations \t" + args[0]);

			long start = System.currentTimeMillis();
			// Initialize population
			Population population = new Population();
			InitializeChromosome initializeChromosome = new InitializeChromosome();
			population.population = initializeChromosome.initiaizePopulation(ReadFiles.clusters, populationLength,
					MyRandom.r);

			for (int i = 0; i < generation; i++) {

				Population subPop = new Population(); // store new
														// individuals are
														// generated
				double[] popFitness = population.populationFitness1();

				for (int j = 0; j < populationLength / 2; j++) {
					Individual father = new Individual();
					Individual mother = new Individual();

					// use tournament to select Individual
					father = population.getIndividual(selection.touranmentSelection(popFitness, 2, MyRandom.r));
					mother = population.getIndividual(selection.touranmentSelection(popFitness, 2, MyRandom.r));

					Individual offspring1 = new Individual();
					Individual offspring2 = new Individual();

					// new random number [0;1]
					double d = MyRandom.r.nextDouble();
					// if random number > crossover rate, add the individual
					// to population without crossover
					if (d < p_rate) {
						// else do crossover, generate two offsprings, then
						// do mutation
						ArrayList<int[]> childs = new ArrayList<int[]>();
						childs = crossover.pruferNumberCrossover(father.getGene1(), mother.getGene1(),
								mfeaPruferNumber.calCulateGenLength(ReadFiles.clusters), MyRandom.r);
						offspring1.setGene(childs.get(0));
						offspring2.setGene(childs.get(1));

					} else {
						offspring1.setGene(mutation.mutationPrufer(ReadFiles.clusters, MyRandom.r, father.getGene1(),
								evaluation.getStartPoint(ReadFiles.clusters), m_rate));
						offspring2.setGene(mutation.mutationPrufer(ReadFiles.clusters, MyRandom.r, mother.getGene1(),
								evaluation.getStartPoint(ReadFiles.clusters), m_rate));

					}
					subPop.addIndiv(offspring1);
					subPop.addIndiv(offspring2);

				}

				// get the index of the best Individual of old Population
				double bestFiness = popFitness[0];
				int bestFinessIndex = 0;
				for (int t = 0; t < populationLength; t++) {
					if (popFitness[t] < bestFiness) {
						bestFiness = popFitness[t];
						bestFinessIndex = t;
					}

				}
				pw3.println(i + "\t" + popFitness[bestFinessIndex]);

				// add the best Individual of old Population to new
				// Population

				subPop.population.set(MyRandom.r.nextInt(populationLength), population.population.get(bestFinessIndex));
				population = subPop;
			}
			pw3.close();
			long end = System.currentTimeMillis();
			double[] popFitness = population.populationFitness1();
			double bestFiness = popFitness[0];
			int bestFinessIndex = 0;
			for (int t = 1; t < populationLength; t++) {
				if (popFitness[t] < bestFiness) {
					bestFiness = popFitness[t];
					bestFinessIndex = t;
				}
			}

			// print to files
			// System.out.println(popFitness[bestFinessIndex]);
			// pw.println(popFitness[bestFinessIndex]);

			System.out.println("|" + loop + "\t|" + args[0] + "\t|" + popFitness[bestFinessIndex] + "\t|"
					+ MEF2Instances.getDateFromMillis(end - start) + "\t|");
			pw.println("Name: " + args[0]);
			pw.println("Seed: " + loop);
			pw.println("Fitness: " + popFitness[bestFinessIndex]);
			pw.println("Time: " + MEF2Instances.getDateFromMillis(end - start));
			// for(int m = 0; m < ReadFiles.num_vertex; m ++){
			//// pw.println();
			//// for( int n = 0; n < ReadFiles.num_vertex; n++){
			//// if(
			// population.getIndividual(bestFinessIndex).getGene()[m][n] >
			// 0){
			//// population.getIndividual(bestFinessIndex).getGene()[m][n] =
			// 1.0f; }
			//// pw.print((int)
			// population.getIndividual(bestFinessIndex).getGene()[m][n]+"\t");
			////
			//// }
			// }
			pw.close();

		}

	}
	public static void sortByCostIndex(ArrayList<Individual> pop, int pop_size, int index) {
		for (int i = 0; i < pop_size; i++) {
			pop.get(i).cost = pop.get(i).getFactorialCost()[index];
		}
		Collections.sort(pop, ChromosomeCmp.compareByFactorialCost);

	}
}

