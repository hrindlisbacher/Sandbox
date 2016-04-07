BinoCR <- function(n,alpha=0.05,
	method=c("exact","wilson","wald"),
	type=c("twosided","upper","lower"),
	display=TRUE)
# BinoCR liefert einen Konfidenzbereich fuer einen
# Binomialparameter p, basierend auf Beob. H ~ Bin(n,p),
# in Form einer Matrix "ab" mit zwei Spalten. Dabei ist
#    ab[j,] = Konfidenzintervall [a(j-1),b(j-1)].
{
	method <- match.arg(method,c("exact","wilson","wald"))
	type <- match.arg(type,c("twosided","upper","lower"))
	ab <- cbind(a=rep(NA,n+1),b=rep(NA,n+1))
	dimnames(ab)[[1]] <- 0:n
	ph <- (0:n)/n
	c1 <- qnorm(1 - alpha)/sqrt(n)
	c2 <- qnorm(1 - alpha/2)/sqrt(n)
	if (method == "exact")
	{
		# Konfidenzintervalle nach Clopper-Pearson:
		for (h in 0:n)
		{
			if (type == "lower")
			{
				ab[h+1,] <- binom.test(h,n,
					conf.level=1-alpha,
					alternative="greater")$conf.int
			}
			if (type == "upper")
			{
				ab[h+1,] <- binom.test(h,n,
					conf.level=1-alpha,
					alternative="less")$conf.int
			}
			if (type == "twosided")
			{
				ab[h+1,] <- binom.test(h,n,
					conf.level=1-alpha,
					alternative="two.sided")$conf.int
			}
		}
	}
	if (method == "wilson")
	{
		# Konfidenzintervalle nach Wilson:
		if (type == "lower")
		{
			ab[,1] <- (ph + c1^2/2 -
				c1*sqrt(ph*(1-ph) + c1^2/4)) /
				(1 + c1^2)
			ab[,2] <- 1
			# Wegen numerischer Fehler:
			ab[1,1] <- 0
		}
		if (type == "upper")
		{
			ab[,1] <- 0
			ab[,2] <- (ph + c1^2/2 +
				c1*sqrt(ph*(1-ph) + c1^2/4)) /
				(1 + c1^2)
			# Wegen numerischer Fehler:
			ab[n+1,2] <- 1
		}
		if (type == "twosided")
		{
			ab[,1] <- (ph + c2^2/2 -
				c2*sqrt(ph*(1-ph) + c2^2/4)) /
				(1 + c2^2)
			ab[,2] <- (ph + c2^2/2 +
				c2*sqrt(ph*(1-ph) + c2^2/4)) /
				(1 + c2^2)
			# Wegen numerischer Fehler:
			ab[1,1] <- 0
			ab[n+1,2] <- 1
		}
	}
	if (method == "wald")
	{
		# Konfidenzintervalle nach Wald:
		if (type == "lower")
		{
			ab[,1] <- pmax(ph - c1*sqrt(ph*(1-ph)), 0)
			ab[,2] <- 1
		}
		if (type == "upper")
		{
			ab[,1] <- 0
			ab[,2] <- pmin(ph + c1*sqrt(ph*(1-ph)), 1)
		}
		if (type == "twosided")
		{
			ab[,1] <- pmax(ph - c2*sqrt(ph*(1-ph)), 0)
			ab[,2] <- pmin(ph + c2*sqrt(ph*(1-ph)), 1)
		}
	}
	if (display)
	{
		# Stelle den Konfidenzbereich graphisch dar:
		plot(0:n,ph,pch=1,xlab="h",ylab="p",main=method)
		for (h in 0:n)
		{
			lines(c(h,h),c(0,1),col="red",lwd=1)
			lines(c(h,h),ab[h+1,],col="darkgreen",lwd=3)
		}
	}
	return(ab)
}


BinoCovProb.naive <- function(ab,alpha=0.05,ylim=c(0,1))
# Fuer einen gegebenen Konfidenzbereich, beschrieben
# durch eine Matrix ab wie in BinoCR(), zeichne seine
# Ueberdeckungswahrscheinlichkeit
#    Pr_p(a(H) <= p <= b(H))
# als Funktion von p in [0,1]. Mit dem optionalen
# Parameter ylim soll die y-Achse allenfalls auf ein
# Teilintervall von [0,1] eingeschraenkt werden. 
{
	n <- dim(ab)[1] - 1
	h <- 0:n
	pp <- sort(c(seq(0,1,by=0.001),ab))
	UebWs <- rep(NA,length(pp))
	for (i in 1:length(pp))
	{
		fh <- dbinom(h,n,pp[i])
		UebWs[i] <- sum(fh[pp[i] >= ab[h+1,1] &
			pp[i] <= ab[h+1,2]])
	}
	plot(c(0,1),rep(1-alpha,2),type="l",col="red",
		ylim=ylim,xlab="p",ylab="Ws(p in [a(H), b(H)])")
	lines(pp,UebWs,lwd="2")
}


BinoCovProb <- function(ab,alpha=0.05,ylim=c(0,1))
# Fuer einen gegebenen Konfidenzbereich, beschrieben
# durch eine Matrix ab wie in BinoCR(), zeichne seine
# Ueberdeckungswahrscheinlichkeit
#    Pr_p(a(H) <= p <= b(H))
# als Funktion von p in [0,1]. Mit dem optionalen
# Parameter ylim soll die y-Achse allenfalls auf ein
# Teilintervall von [0,1] eingeschraenkt werden. 
{
	# Bestimme Stichprobenumfang:
	n <- dim(ab)[1] - 1
	# Moegliche Werte von H:
	h <- 0:n
	
	# Feines Gitter von Parametern p, fuer welche
	# die Ueberdeckungswahrscheinlichkeit berechnet
	# werden soll:
	pp <- sort(c(seq(0,1,by=0.001),ab))
	# Entferne doppelte Eintraege:
	pp <- unique(pp)
	
	# Ueberdeckungswahrscheinlichkeiten sowie
	# links- und rechtsseitige Grenzwerte:
	UebWs <- rep(NA,length(pp))
	UebWs.l <- rep(NA,length(pp))
	UebWs.r <- rep(NA,length(pp))
	for (i in 1:length(pp))
	{
		fh <- dbinom(h,n,pp[i])
		UebWs[i]   <- sum(fh[pp[i] >= ab[h+1,1] &
		                     pp[i] <= ab[h+1,2]])
		UebWs.l[i] <- sum(fh[pp[i] >  ab[h+1,1] &
		                     pp[i] <= ab[h+1,2]])
		UebWs.r[i] <- sum(fh[pp[i] >= ab[h+1,1] &
		                     pp[i] <  ab[h+1,2]])
	}
	# Nun verdreifache jede Komponente von pp, damit
	# man jeweils den links-seitigen Grenzwert, den Wert
	# und den rechtsseitigen Grenzwert der Uebrderckungs-
	# wahrscheinlichkeit angeben kann:
	PP <- rep(pp,each=3)
	UEBWS <- as.vector(rbind(UebWs.l,UebWs,UebWs.r))
	# Entferne den linksseitigen Grenzwert bei 0 und
	# den rechtsseitigen Grenzwert bei 1:
	PP <- PP[2:(length(PP)-1)]
	UEBWS <- UEBWS[2:(length(UEBWS)-1)]
	plot(c(0,1),rep(1-alpha,2),type="l",col="red",
		ylim=ylim,xlab="p",ylab="Ws(p in [a(H), b(H)])")
	lines(PP,UEBWS,lwd="2")
}
