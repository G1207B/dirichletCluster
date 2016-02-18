package edu.ucr.cuilab.algorithm;

import java.util.Comparator;

public class DoubleReadCompareGC implements Comparator<DoubleRead> {
	@Override
	public int compare(DoubleRead arg0, DoubleRead arg1) {
		if (arg0.getGC() > arg1.getGC()) {
			return 1;
		} else if (arg0.getGC() < arg1.getGC()){
			return -1;
		} else {
			return 0;
		}
	}
}
