package operator;

import java.util.Comparator;
import structures.Cluster;
import structures.Cycles;
import structures.Individual;

public class ChromosomeCmp implements Comparator<Individual> {

	@Override
	public int compare(Individual ind1, Individual ind2) {
		// TODO Auto-generated method stub
		return (ind1.getScalarFitness() < ind2.getScalarFitness() ? -1
				: ind1.getScalarFitness() < ind2.getScalarFitness() ? 1 : 0);
	}

	public static Comparator<Individual> compareByFactorialCost = new Comparator<Individual>() {
		public int compare(Individual other, Individual one) {
			return Double.compare(other.cost, one.cost);
		}
	};
	public static Comparator<Individual> compareByScalarFitness = new Comparator<Individual>() {
		public int compare(Individual one, Individual other) {
			return Double.compare(other.getScalarFitness(), one.getScalarFitness());
		}
	};
	public static Comparator<Cycles> compareByMaxElement = new Comparator<Cycles>() {
		public int compare(Cycles cc1, Cycles cc2) {
			return Integer.compare(cc2.getMaxElement(), cc1.getMaxElement());
		}
	};
	public static Comparator<Cluster> compareByNumberOfCluster = new Comparator<Cluster>() {
		public int compare(Cluster cluster1, Cluster cluster2) {
			return Integer.compare(cluster2.getCluster().size(), cluster1.getCluster().size());
		}
	};

}
