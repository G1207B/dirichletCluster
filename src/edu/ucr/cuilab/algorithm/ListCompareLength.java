package edu.ucr.cuilab.algorithm;

import java.util.Comparator;
import java.util.Set;

public class ListCompareLength implements Comparator<Set<Integer>> {

	@Override
	public int compare(Set<Integer> arg0, Set<Integer> arg1) {
		if (arg0.size() > arg1.size()) {
			return -1;
		} else {
			return 1;
		}
	}

}
