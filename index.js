// Coeficientes baseados nas variáveis aleatórias
var A, B, C, D, E, F, G, H;

// Valores iniciais da equação diferencial
var ini_y , ini_y_d;

// Funções vetoriais
// Multiplicar por escalar
var scale = function(scalar, vector){
	var v = [];
	for ( var i = 0; i < vector.length; i++ ){
		v.push(scalar * vector[i]);
	}
	return v;
};
// Produto interno
var innerProduct = function(v1, v2){
	var p = 0;
	for ( var i = 0; i < v1.length; i++ ){
		p += v1[i] * v2[i];
	}
	return p;
};
// Tamanho do vetor
var norm = function(vector){
	return Math.sqrt(innerProduct(vector, vector));
};

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

// Variável aleatória exponencial
var expon = function(lambda){
	return Math.log(1-Math.random())/(-lambda);
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

// Densidade de probabilidade f_beta (para a espiral)
// x: Valor no intervalo [0, 10]
var f_beta = function(x){
	var c = 3.60566;
	return c/(Math.exp(2*x)*Math.pow((1+x), 2));
}

// Ponto uniformemente distribuído na espiral
var Espiral = function(){
	// Gerar variável aleatória com densidade 
	// f_beta pelo método aceitação/rejeição
	var x, y;
	x = 10*Math.random();
	y = 4*Math.random();
	while ( f_beta(x) < y ){
		x = 10*Math.random();
		y = 4*Math.random();
	}
	var beta = x;
	// Ângulo uniformemente distribuído entre -2*pi e Beta
	var teta = Math.log((Math.exp(beta)-Math.exp(-2*Math.PI))*Math.random() + Math.exp(-2*Math.PI));
	// Distância da origem para para o ângulo alfa
	var rho = Math.exp(teta);
	// Coeficiente E
	E = -rho;
	// Coeficiente G
	G = teta;
};

// Ponto uniformemente distribuído no sólido (funil)
var Funil = function(){
	// Ponto dentro do disco na altura z 
	var z = 0.5*Math.log(1-Math.random())
	// Maior raio do disco que corta o funil
	var rz = Math.exp(z);
	// Gerar coordenadas x,y dentro do disco
	var x = -rz + 2*Math.random()*rz;
	var y = -rz + 2*Math.random()*rz;
	// Gerar novamente caso ponto fora do disco ou sem tamanho
	while ( Math.pow(x,2) + Math.pow(y,2) > Math.pow(rz,2) || ( x==0 && y==0 ) ){
		x = -rz + 2*Math.random()*rz;
		y = -rz + 2*Math.random()*rz;
	}
	var teta = Math.atan2(y, x);
	if ( teta < 0 ){
		// Ângulo negativo implica que está entre pi e 2*pi
		// Corrigir para continuar no sentido anti-horário
		teta = 2*Math.PI + teta;
	}
	// Valor esperado para o raio do disco
	// calculado numericamente pela média em
	// 1.000.000.000 de amostras
	var ER_v = 4/9;
	// Distância do ponto ao eixo z
	var R_v = Math.sqrt(Math.pow(x, 2) + Math.pow(y, 2));
	// Coeficiente C
	C = 1 + R_v - ER_v;
	// Coeficiente D
	D = teta;
	// Coeficiente F
	F = z;
	// Condição inicial y(0)
	ini_y = (teta - Math.PI)/Math.PI;
};

// Ponto uniformemente distribuído na superfície da esfera
var Esfera = function(){
	// Raio da esfera
	var rho = 1;
	// Armazenar par normal gerado pelo método Box-Muller
	var par = [0,0];
	// Posicao
	var ponto = [0,0,0];
	par = parNormal();
	// x
	ponto[0] = par[0];
	// y
	ponto[1] = par[1];
	par = parNormal();
	// z
	ponto[2] = par[0];
	// Projetar na superfície da esfera
	var t = norm(ponto);
	ponto = scale(rho/t, ponto);
	// Ângulo teta: 0 < teta < 2*pi
	var teta = Math.atan2(ponto[1], ponto[0]);
	if ( teta < 0 ){
		// Ângulo negativo implica que está entre pi e 2*pi
		// Corrigir para continuar no sentido anti-horário
		teta = 2*Math.PI + teta;
	}
	// Ângulo phi: -pi/2 < phi < pi/2
	// Função Math.acos retorna valor [0,pi], necessário rotacionar
	var phi = Math.acos(ponto[2]) - Math.PI/2;
	// Coeficiente A
	A = 1 + teta - Math.PI;
	// Coeficiente B
	B = 1 + phi;
	// Condição inicial y´(0)
	ini_y_d = (teta - Math.PI)/Math.PI;
};

// Gerar coeficientes a partir dos pontos aleatoriamente distribuídos
// exibir: Exibir valores nos campos das páginas (0 ou 1)
var coeficientes = function(exibir){
	Espiral();
	Funil();
	Esfera();
	// Coeficiente H
	// "Grudar" todos os intervalos de C (a união) e identifica-los com sub-intervalos de (0, 1).
	var x = Math.random();
	var p = 1/2;
	var i = 0;
	// Iterar até identificar o intervalo de C
	while ( p < x ){
		i++;
		p += 1/Math.pow(2, i+1);
	}
	// Mapear H de acordo com o intervalo onde x caiu
	p = p - 1/Math.pow(2, i+1);
	H = i + x - p;

	if (exibir == 1){
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
	}
};

// Parcela M(t) da solução não-homogênea
var M = function(t){
	// Parâmetro do processo de Poisson
	var lambda_p = geom(1/5);
	// Gerar lambda_p uniformes em (0,1)
	var S = [];
	for (var s = 0; s < lambda_p; s++){
		S.push(Math.random());
	}
	// Ordenar eventos
	S.sort();
	// Contar quantos eventos do processo de poisson ocorreram até t
	var c = 0;
	for (var s = 0; s < lambda_p; s++){
		if ( S[s] < t ){
			// Somar evento
			c++;
		}
		else {
			// Instante do evento ultrapassou t
			break;
		}
	}
	// Somatória
	var sum = 0;
	for (var i = 0; i < c; i++){
		// Variável z = 1 com probabilidade 1/2 e 0 com prob. 1/2
		var z = ( Math.random() < .5 ) ? 1 : 0;
		var x = expon(1+poiss(.5)*(1+t));
		sum += Math.pow(-1, z)*x;
	}
	return sum;
}

var ler_x = function(){
	var x = document.querySelector("input[name='x']").value;
	// Verificar se x é válido
	if (isNaN(x)) x = 0;
	// Não permitir valores fora do intervalo [0,1]
	if ( x < 0 ) x = 0;
	if ( x > 1 ) x = 1;
	document.getElementById("x-in-sum").textContent = x;
	return x;
}

// Gerar solução para a equação não-homogênea a partir 
// dos coeficientes e condições iniciais aleatórios.
// x: Variável x em [0,1]
// exibir: Mostrar valores nos campos da página (0 ou 1)
var solucao = function(x, exibir){
	// Valor da somatória (solução particular iii) para o x escolhido
	var m = M(x);
	if (exibir == 1){
		document.getElementById("sum").textContent = m;
	}

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

	// Segunda solução particular (seno e cosseno)
	// Cosseno
	sol_p[1] = ((F*(C-A*Math.pow(G,2)))/(Math.pow(C-A*Math.pow(G,2),2)+Math.pow(2*B*G,2)))*Math.cos(G*x+H);
	// Seno
	sol_p[1] += ((2*B*G*F)/(Math.pow(C-A*Math.pow(G,2),2)+Math.pow(2*B*G,2)))*Math.sin(G*x+H);

	// Terceira solução particular (somatória)
	sol_p[2] = m/C;

	// Resultado
	if (exibir == 1){
		document.getElementById("y").textContent = c[0]*sol[0] + c[1]*sol[1] + sol_p[0] + sol_p[1] + sol_p[2];
	}
	return c[0]*sol[0] + c[1]*sol[1] + sol_p[0] + sol_p[1] + sol_p[2];
}

// Computar solução
document.querySelector("button[name='solve']").addEventListener("click", function(event){
	solucao(ler_x(), 1);
});

// Gerar coeficientes a partir dos pontos aleatoriamente distribuídos
document.querySelector("button[name='coef']").addEventListener("click", function(event){
	coeficientes(1);
	// Calcular a solução
	solucao(ler_x(), 1);
});

// Gerar coeficientes e computar solução ao carregar a página
coeficientes(1);
// Calcular a solução
solucao(ler_x(), 1);


/*
 * Funções para estimar esperança e variância de y(1)
 * são invocadas maualmente a partir do console
 */

// Calcular esperança para y(1)
// n: Número de iterações
var esperanca = function(n){
	// Soma dos resultados
	var sum = 0;
	// Cópia da solução
	var sol = 0;
	// Cópia do número de iterações
	var d = n;
	for (var i = 0; i < n; i++){
		// Gerar coeficientes sem exibir na página
		coeficientes(0);
		// Calcular o resultado sem exibir na página
		sol = solucao(1, 0);
		if ( sol == Number.POSITIVE_INFINITY || sol == Number.NEGATIVE_INFINITY || isNaN(sol) ){
			// Erro na solução numérica, descartar essa iteração
			d--;
			continue;
		}
		sum += sol;
	}
	// Média para y(1)
	console.log(d, sum / d);

}

// Calcular variância para y(1)
// e: esperança
// n: Número de iterações
var variancia = function(e, n){
	// Soma dos resultados
	var sum = 0;
	// Cópia da solução
	var sol = 0;
	// Cópia do número de iterações
	var d = n;
	for (var i = 0; i < n; i++){
		// Gerar coeficientes sem exibir na página
		coeficientes(0);
		// Calcular o resultado sem exibir na página
		sol = solucao(1, 0);
		if ( sol == Number.POSITIVE_INFINITY || sol == Number.NEGATIVE_INFINITY || isNaN(sol) ){
			// Erro na solução numérica, descartar essa iteração
			d--;
			continue;
		}
		// Variância
		sum += Math.pow(sol - e, 2);
	}
	// Média para y(1)
	console.log(d, sum / d);

}
