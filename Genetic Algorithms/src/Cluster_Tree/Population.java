package Cluster_Tree;

import java.util.ArrayList;
import java.util.Random;

public class Population {
	public static int populationLength = 100; // default length of population
	public ArrayList<Individual> population = new ArrayList<Individual>(populationLength);

	public Population() {
	};
	
	// get gene of an Individual
	public double[][] getIndividual(int i) {
		return population.get(i).getGene();
	}
	

	
	public void addIndiv(Individual individual){
		population.add(individual);
	}

	
	}
