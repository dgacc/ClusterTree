package Cluster_Tree;

import java.util.Random;

public class Selection {
	Population  population = new Population();
	public Individual selectIndiv(int tournament_size, Random r) {
		int k = r.nextInt(Population.populationLength);
		for (int i = 0; i < tournament_size; i++) {
			int j = r.nextInt(Population.populationLength);
			if (population.getIndividual(k).getClass(). > population.get(j).getFitness()) {
				k = j;
			}
		}
		return population.get(k);
	}


}
