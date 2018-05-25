// Coeficientes baseados nas variáveis aleatórias
var A, B, C, D, E, F, G, H;

// Valores iniciais da equação diferencial
var ini_y , ini_y_d;

// Variável aleatória geométrica (Número de falhas até o primeiro sucesso)
var geom = function(lambda){
	return Math.floor(Math.log(1-Math.random())/Math.log(1-lambda));
}

// Variável aleatório Poisson
var poiss = function(lambda){
	var L = Math.exp(-lambda);
	var k = 0;
	var p = 1;
	k++;
	p *= Math.random();
	while ( p > L ){
		k++;
		p *= Math.random();
	}
	return k - 1;
}

// Um par de variáveis aleatórias com distribuição normal,
// esperança zero e variância unitária (método Box-Muller)
var parNormal = function(){
	// Gerar um ponto uniforme dentro do disco 
	// de raio unitário centrado na origem
	var p = [-1 + 2*Math.random(), -1 + 2*Math.random()];
	var r = Math.pow(p[0],2) + Math.pow(p[1],2);
	// Gerar novamente caso ponto fora 
	// do disco ou exatamente na origem
	while ( r > 1 || r == 0 ){
		p = [-1 + 2*Math.random(), -1 + 2*Math.random()];
		r = Math.pow(p[0],2) + Math.pow(p[1],2);
	}
	// Fator comum
	var f = Math.sqrt(-2*Math.log(r)/r);
	// Par normal
	return [p[0]*f,p[1]*f];
}

// Pontos uniformemente distribuídos na espiral
var uniformeEspiral = function(){
	// Ponto uniformemente distribuído na espiral entre 0 e 2*PI.
	// Não incluir 2*PI para não gerar raio de tamanho zero
	var beta = 2*Math.PI - 1e-10;
	var a = 1 - Math.random();
	// Ângulo do ponto
	var alfa = beta + Math.log(a);
	// Raio do ponto
	var rho = Math.exp(beta) * a;
	// Coordenadas cartesianas
	var x = rho * Math.cos(alfa);
	var y = rho * Math.sin(alfa);
	// Posicao
	var ponto = [x,y,0];
	// Coeficiente E
	E = -rho;
	// Coeficiente G
	G = alfa;
	// Coeficiente H
	H = Math.sqrt(2)*(Math.exp(alfa)-1);
};

// Pontos uniformemente distribuídos no sólido
var uniformeSolido = function(){
	// Ponto dentro do disco na altura z 
	var z = 0.5*Math.log(1-Math.random())
	// Raio do disco
	var ez = Math.exp(z);
	var x = -ez + 2*Math.random()*ez;
	var y = -ez + 2*Math.random()*ez;
	// Gerar novamente caso ponto fora do disco ou sem tamanho
	while ( Math.pow(x,2) + Math.pow(y,2) > Math.pow(ez,2) || ( x==0 && y==0 ) ){
		x = -ez + 2*Math.random()*ez;
		y = -ez + 2*Math.random()*ez;
	}
	// Projetar na circunferência
	var t = Math.sqrt(Math.pow(x,2)+Math.pow(y,2));
	x *= ez/t;
	y *= ez/t;
	// Posicao
	var ponto = [x,y,z];
	// Valor esperado para o raio do disco
	var esp = 0.5;
	// Coeficiente C
	C = 1 + ez - esp;
	// Coeficiente D
	D = Math.atan2(ponto[1], ponto[0]);
	// Coeficiente F
	F = z;
	// Condição inicial y(0)
	ini_y = (Math.atan2(ponto[1], ponto[0]) - Math.PI)/Math.PI;
};

// Pontos uniformemente distribuídos na superfície da esfera
var uniformeEsfera = function(){
	// Raio da esfera
	var rho = 1;
	// Armazenar par normal gerado pelo método Box-Muller
	var par = [0,0];
	// Posicao
	var ponto = [0,0,0];
	par = parNormal();
	ponto[0] = par[0];
	ponto[1] = par[1];
	par = parNormal();
	ponto[2] = par[0];
	// Projetar na superfície da esfera
	var t = Op.norm(ponto);
	ponto = Op.scale(rho/t, ponto);
	// Coeficiente A
	A = 1 + 2*Math.acos(ponto[2]) - Math.PI;
	// Coeficiente B
	B = Math.atan(ponto[1]/ponto[0]) + 1;
	// Condição inicial y´(0)
	ini_y_d = (2*Math.acos(ponto[2]) - Math.PI)/Math.PI;
};

// Gerar coeficientes a partir dos pontos aleatoriamente distribuídos
var coeficientes = function(){
	uniformeEspiral();
	uniformeSolido();
	uniformeEsfera();
	document.querySelector("span[data-variable='A']").textContent = A;
	document.querySelector("span[data-variable='B']").textContent = B;
	document.querySelector("span[data-variable='C']").textContent = C;
	document.querySelector("span[data-variable='D']").textContent = D;
	document.querySelector("span[data-variable='E']").textContent = E;
	document.querySelector("span[data-variable='F']").textContent = F;
	document.querySelector("span[data-variable='G']").textContent = G;
	document.querySelector("span[data-variable='H']").textContent = H;
	document.querySelector("span[data-ini='y']").textContent = ini_y;
	document.querySelector("span[data-ini='y_d']").textContent = ini_y_d;
};

// Parcela M(t) da solução não-homogênea
var M = function(t){
	// Parâmetro do processo de Poisson
	var lambda_p = geom(1/5);
	// Número de ocorrências do processo de Poisson no intervalo t
	N = poiss(t*lambda_p);
	// Poisson da variável x
	var lambda_e = poiss(.5);
	// Somatória
	var sum = 0;
	for (var i = 0; i < N; i++){
		// Variável z = 1 com probabilidade 1/2 e 0 com prob. 1/2
		var z = ( Math.random() < .5 ) ? 1 : 0;
		var x = Math.exp(1+lambda_e*(1+t));
		sum += Math.pow(-1, z)*x;
	}
	return sum;
}

// Gerar solução para a equação não-homogênea a partir 
// dos coeficientes e condições iniciais aleatórios.
var solucao = function(){
	// Variável x em [0,1]
	var x = document.querySelector("input[name='x']").value;
	// Varificar se x é válido
	if (isNaN(x)) x = 0;
	// Não permitir valores fora do intervalo [0,1]
	if ( x < 0 ) x = 0;
	if ( x > 1 ) x = 1;
	document.getElementById("x-in-sum").textContent = x;

	// Valor da somatória (solução particular iii) para o x escolhido
	var m = M(x);
	document.getElementById("sum").textContent = m;

	// Raízes das soluções elementares
	var r = [0,0];
	// Valor das soluções elementares
	var sol = [0,0];
	// Valor das derivadas das soluções elementares
	var sol_d = [0,0];
	// Constantes da solução homogênea
	var c = [0,0];
	// Evitar divisão por zero
	while ( A == 0 ){
		// Gerar coeficientes aleatórios até obter um A 
		// diferente de zero (provavelmente A nunca será zero)
		coeficientes();
	}
	// Determinar solução homogênea
	var delta = Math.pow(2*B, 2) - 4*A*C;
	if ( delta > 0 ){
		// Soluções elementares para delta positivo
		r[0] = (-2*B+Math.sqrt(delta))/(2*A);
		r[1] = (-2*B-Math.sqrt(delta))/(2*A);
		sol[0] = Math.exp(r[0]*x);
		sol[1] = Math.exp(r[1]*x);
		sol_d[0] = r[0]*Math.exp(r[0]*x);
		sol_d[1] = r[1]*Math.exp(r[1]*x);
	}
	if ( delta == 0 ){
		// Soluções elementares para delta igual a zero
		r[0] = -B/A;
		sol[0] = Math.exp(r[0]*x);
		sol[1] = x*Math.exp(r[0]*x);
		sol_d[0] = r[0]*Math.exp(r[0]*x);
		sol_d[1] = r[0]*x*Math.exp(r[0]*x);
	}
	if ( delta < 0 ){
		// Soluções elementares para delta negativo
		// Usar o vetor de raízes para armazenar os valores
		r[0] = -B/A;
		r[1] = Math.sqrt(Math.abs(delta))/(2*A);
		sol[0] = Math.exp(r[0]*x)*Math.cos(r[1]*x);
		sol[1] = Math.exp(r[0]*x)*Math.sin(r[1]*x);
		sol_d[0] = r[0]*Math.exp(r[0]*x)*Math.cos(r[1]*x) - Math.exp(r[0]*x)*r[1]*Math.sin(r[1]*x);
		sol_d[1] = r[0]*Math.exp(r[0]*x)*Math.sin(r[1]*x) + Math.exp(r[0]*x)*r[1]*Math.cos(r[1]*x);
	}
	// Determinar constantes da solução homogênea a partir dos valores iniciais
	// Solução do sistema 2x2 pela regra de Cramer
	var det = sol[0]*sol_d[1] - sol[1]*sol_d[0];
	var det_c0 = ini_y*sol_d[1] - sol[1]*ini_y_d;
	var det_c1 = sol[0]*ini_y_d - ini_y*sol_d[0];
	c[0] = det_c0/det;
	c[1] = det_c1/det;

	// Soluções particulares
	var sol_p = [0,0,0];

	// Primeira solução particular (exponencial)
	sol_p[0] = (D/(A*Math.pow(E, 2)+2*B*E+C))*Math.exp(E*x);

	// Segunda solução particular (cosseno)
	sol_p[1] = (F/(-A*Math.pow(G, 2)+2*B*G+C))*Math.cos(G*x + H);

	// Terceira solução particular (somatória)
	sol_p[2] = C*m;

	// Resultado
	document.getElementById("y").textContent = c[0]*sol[0] + c[1]*sol[1] + sol_p[0] + sol_p[1] + sol_p[2];
}

// Computar solução
document.querySelector("button[name='solve']").addEventListener("click", function(event){
	solucao();
});

// Gerar coeficientes a partir dos pontos aleatoriamente distribuídos
document.querySelector("button[name='coef']").addEventListener("click", function(event){
	coeficientes();
});

// Gerar coeficientes e computar solução ao carregar a página
coeficientes();
solucao();
