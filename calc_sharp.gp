mzv_basis_uncalculated(d) = {
	if(d<=2, return([
		[[]],
		[],
		[2]
	][d+1]));
	concat(
		apply(s -> concat(3,s), mzv_basis_uncalculated(d-3)),
		apply(s -> concat(2,s), mzv_basis_uncalculated(d-2))
	)
}
mzv_basis(d) = {
	return(apply(s -> zetamult(s), mzv_basis_uncalculated(d)))
}
zeta_star(k) = {
	s = 0;
	i = 0;
	for(i=0,2^(#k-1)-1,
		l = digits(i,2);
		l = concat(vector(#k-1-#l,i,0),l);
		t = vector(vecsum(l)+1, i, 0);
		t_i = 1;
		for(j=1,#k-1,
			t[t_i] += abs(k[j]) * if(t[t_i]<0,-1,1);
			t[t_i] *= sign(k[j]);
			t_i += l[j];
		);
		t[t_i] += abs(k[#k]) * if(t[t_i]<0,-1,1);
		t[t_i] *= sign(k[#k]);
		s += polylogmult(Vec(abs(t)),Vec(apply(sign,t)))
		
	);
	return(s)
}
zeta_sharp(k) = {
	s = 0;
	i = 0;
	for(i=0,2^(#k-1)-1,
		l = digits(i,2);
		l = concat(vector(#k-1-#l,i,0),l);
		t = vector(vecsum(l)+1, i, 0);
		t_i = 1;
		for(j=1,#k-1,
			t[t_i] += abs(k[j]) * if(t[t_i]<0,-1,1);
			t[t_i] *= sign(k[j]);
			t_i += l[j];
		);
		t[t_i] += abs(k[#k]) * if(t[t_i]<0,-1,1);
		t[t_i] *= sign(k[#k]);
		s += 2 ^ (vecsum(l) + 1) * polylogmult(Vec(abs(t)),Vec(apply(sign,t)))
		
	);
	return(s)
}
ramified(k) = vecmax(abs(lindep(concat(zeta_sharp(k),mzv_basis(vecsum(abs(k)))))))

admissible(d) = {
	if(d<=2, return([
		[[]],
		[[-1]],
		[[2],[-2]]
	][d+1]));
	concat(
	concat(
		apply(r -> apply(s -> concat(s,r), admissible(d-r)), [1..d]),
	),
	concat(
		apply(r -> apply(s -> concat(s,-r), admissible(d-r)), [1..d])
	)
	)
}
search(d) = Vecrev(vecsort(select(p -> p[1] < 10 ^ (2 * d), apply(r -> [ramified(r),r], admissible(d))),1))
search_two(d) = Vecrev(vecsort(select(p -> p[1] < 10 ^ (2 * d), apply(r -> [ramified(r),r], apply(k->concat(k,2),admissible(d-2)))),1))
