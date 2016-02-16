package edu.ucr.cuilab.algorithm;

import java.util.Comparator;
import java.util.List;

public class ListCompareLength implements Comparator<List<Integer>> {

	@Override
	public int compare(List<Integer> arg0, List<Integer> arg1) {
		if (arg0.size() > arg1.size()) {
			return -1;
		} else {
			return 1;
		}
	}

}
