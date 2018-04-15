package operator;
import java.util.Comparator;

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
}
