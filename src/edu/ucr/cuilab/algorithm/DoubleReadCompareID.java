package edu.ucr.cuilab.algorithm;

import java.util.Comparator;

public class DoubleReadCompareID implements Comparator<DoubleRead> {
	@Override
	public int compare(DoubleRead arg0, DoubleRead arg1) {
		if (arg0.getId() > arg1.getId()) {
			return 1;
		} else {
			return -1;
		}
	}
}
