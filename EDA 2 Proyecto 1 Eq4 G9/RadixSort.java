import java.io.FileWriter;
import java.io.FileReader;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.File;
import java.util.Scanner;

/** La clase RadixSort se encarga de realizar el ordenamiento externo de un archivo bajo un formato específico.
    Esta clase unicamente realiza el ordenamiento por los numero de claves. Para su funcionamiento se utilizan diferentes
    utilerias de Java. Esta inspirado en el algoritmo Radix Sort interno, solamente extrapolado para ordenamiento externo.
    @author David Calderon
    @author Humberto Hernandez
    @author Inaky Ordiales

*/
public class RadixSort{
	
	private File archivoDesordenado;
	private BufferedReader reader;


	/**
		Metodo para elaborar el ordenamiento externo. En el se checa primero si el nombre del archivo que se le paso
		lleva a un archivo existente. En caso afirmativo crea los archivos auxiliares de colas (metodo crearArchivosColas) y 
		llama al metodo radix, donde se lleva a cabo todo el proceso de analisis y queda el archivo ordenado.
		@param nombreArchivo Es la direccion del archivo a ordenar.
		@return true si se trabaja correctamente el archivo, false si no se encuentra el archivo.
	*/
	public boolean sort(String nombreArchivo){

		if (setArchivoDesordenado(nombreArchivo)){
			crearArchivosColas();
			radix(nombreArchivo);
			return true;
		}else{
			return false;
		}
	}

	/**
		El metodo radix realiza el ordenamiento para radix sort. Se realiza un ciclo de 6 vueltas porque el numero 
		de digitos de las cuentas es 6, y en cada vuelta se lee linea por linea el archivo. Cada linea se pasa para 
		que se clasifique y escriva en el archivo cola auxiliar correspondiente y luego se vuelve a armar el archivo
		original al vaciar cada cola en orden de la 0 a la 9.
		@param nombreArchivo Es la direccion del archivo a ordenar.
	*/
	public void radix(String nombreArchivo){
		
		int k;
		String iter = "Iteraciones_Radix/Iter_Radix_"+nombreArchivo;
		File directorio = new File("Iteraciones_Radix");
		directorio.mkdir();
		File itera = new File("Iteraciones_Radix\\Iter_Radix_"+nombreArchivo);

		//Se borra el contenido del archivo iter en caso de que ya se haya creado en otra ejecucion.
		try{
			FileWriter iteraciones = new FileWriter(iter);
		}catch (IOException e){
			e.printStackTrace();
		}

		for (k = 0; k<6; k++){
			try{
				FileWriter iteraciones = new FileWriter(iter, true);
				reader = new BufferedReader(new FileReader(nombreArchivo));
				String line = reader.readLine();
				while(line != null && line != "\n"){
					clasificarMiembro(line, k);
					line = reader.readLine();
				}
				reader.close();
				iteraciones.append("\n\n\n\n**** ITERACION "+(k+1)+":\n\n");
				regresarArchivo(nombreArchivo, iteraciones);
				iteraciones.close();
				limpiarColas();
			}catch(IOException e){
				e.printStackTrace();
			}
		}

		
	}

	/** 
		Este metodo recibe la linea del archivo a ordenar (con el nombre, apellido, #cuenta) y el entero k, que nos dira el 
		numero de digito de derecha a izquierda que se clasificara.
		@param line Es una linea de un archivo donde se encuentra el numero de cuenta, el cual se analizara.
		@param k El digito del numero de cuenta que se clasificara.
	*/
	public void clasificarMiembro(String line, int k){
		
		int size = line.length();
		int i = 0;
		FileWriter escritura = null;

		String datos[] = line.split(",");

		String nombre = datos[0];
		String apellidos = datos[1];
		String cuenta = datos[2];

		char digito = cuenta.charAt(cuenta.length()-1-k);

		try{
			switch(digito){
				case '0':
					escritura = new FileWriter("Temporal/queue0.txt", true);
					escritura.append(line+"\n");
					escritura.close();
					break;
				case '1':
					escritura = new FileWriter("Temporal/queue1.txt", true);
					escritura.append(line+"\n");
					escritura.close();
					break;
				case '2':
					escritura = new FileWriter("Temporal/queue2.txt", true);
					escritura.append(line+"\n");
					escritura.close();
					break;
				case '3':
					escritura = new FileWriter("Temporal/queue3.txt", true);
					escritura.append(line+"\n");
					escritura.close();
					break;
				case '4':
					escritura = new FileWriter("Temporal/queue4.txt", true);
					escritura.append(line+"\n");
					escritura.close();
					break;
				case '5':
					escritura = new FileWriter("Temporal/queue5.txt", true);
					escritura.append(line+"\n");
					escritura.close();
					break;
				case '6':
					escritura = new FileWriter("Temporal/queue6.txt", true);
					escritura.append(line+"\n");
					escritura.close();
					break;
				case '7':
					escritura = new FileWriter("Temporal/queue7.txt", true);
					escritura.append(line+"\n");
					escritura.close();
					break;
				case '8':
					escritura = new FileWriter("Temporal/queue8.txt", true);
					escritura.append(line+"\n");
					escritura.close();
					break;
				case '9':
					escritura = new FileWriter("Temporal/queue9.txt", true);
					escritura.append(line+"\n");
					escritura.close();
					break;
			}
		}catch(IOException e){
				e.printStackTrace();
		}
	}

	/**
		Vacia todas las colas auxiliares en el archivo original.
		@param nombreArchivo El nombre del archivo que se ordena.
		@param iteraciones El archivo donde se imprime el resultado de cada iteracion.
	*/
	public void regresarArchivo(String nombreArchivo, FileWriter iteraciones){

		try{

			FileWriter escritura = new FileWriter(nombreArchivo);

			reader = new BufferedReader(new FileReader("Temporal/queue0.txt"));
			String line = reader.readLine();
			while(line != null && line != "\n"){
				escritura.append(line+"\n");
				iteraciones.append(line+"\n");
				line = reader.readLine();
			}
			reader.close();
			reader = new BufferedReader(new FileReader("Temporal/queue1.txt"));
			line = reader.readLine();
			while(line != null && line != "\n"){
				escritura.append(line+"\n");
				iteraciones.append(line+"\n");
				line = reader.readLine();
			}
			reader.close();
			reader = new BufferedReader(new FileReader("Temporal/queue2.txt"));
			line = reader.readLine();
			while(line != null && line != "\n"){
				escritura.append(line+"\n");
				iteraciones.append(line+"\n");
				line = reader.readLine();
			}
			reader.close();
			reader = new BufferedReader(new FileReader("Temporal/queue3.txt"));
			line = reader.readLine();
			while(line != null && line != "\n"){
				escritura.append(line+"\n");
				iteraciones.append(line+"\n");
				line = reader.readLine();
			}
			reader.close();
			reader = new BufferedReader(new FileReader("Temporal/queue4.txt"));
			line = reader.readLine();
			while(line != null && line != "\n"){
				escritura.append(line+"\n");
				iteraciones.append(line+"\n");
				line = reader.readLine();
			}
			reader.close();
			reader = new BufferedReader(new FileReader("Temporal/queue5.txt"));
			line = reader.readLine();
			while(line != null && line != "\n"){
				escritura.append(line+"\n");
				iteraciones.append(line+"\n");
				line = reader.readLine();
			}
			reader.close();
			reader = new BufferedReader(new FileReader("Temporal/queue6.txt"));
			line = reader.readLine();
			while(line != null && line != "\n"){
				escritura.append(line+"\n");
				iteraciones.append(line+"\n");
				line = reader.readLine();
			}
			reader.close();
			reader = new BufferedReader(new FileReader("Temporal/queue7.txt"));
			line = reader.readLine();
			while(line != null && line != "\n"){
				escritura.append(line+"\n");
				iteraciones.append(line+"\n");
				line = reader.readLine();
			}
			reader.close();
			reader = new BufferedReader(new FileReader("Temporal/queue8.txt"));
			line = reader.readLine();
			while(line != null && line != "\n"){
				escritura.append(line+"\n");
				iteraciones.append(line+"\n");
				line = reader.readLine();
			}
			reader.close();
			reader = new BufferedReader(new FileReader("Temporal/queue9.txt"));
			line = reader.readLine();
			while(line != null && line != "\n"){
				escritura.append(line+"\n");
				iteraciones.append(line+"\n");
				line = reader.readLine();
			}
			reader.close();
			escritura.close();
			iteraciones.close();
		}catch(IOException e){
				e.printStackTrace();
		}
	}

	/**
		Verifica la existencia de un archivo con el nombre/direccion que se le pasa. En caso de encontrarlo regresa true, 
		si no lo encuentra regresa false.
		@param nombreArchivo El nombre del archivo que se ordena.
		@return true si el archivo ya existe, false si el archivo no existe o no lo encuentra.

	*/
	private boolean setArchivoDesordenado(String nombreArchivo){
		archivoDesordenado = new File(nombreArchivo);
		if (archivoDesordenado.exists()){
			System.out.println("\n El archivo se encontro con exito.");
			return true;
		}else{
			System.out.println("\n No se encuentra el archivo !!!");
			return false;
		}
	}

	/**
		Se crea una carpeta temporal (llamada Temporal) en el lugar de ejecucion del codigo. Luego se crean los 10 archivos 
		de colas auxiliares necesarios para el ordenamiento externo por Radix. Tanto la carpeta como los archivos se les da el
		metodo de que se eliminen al terminar la ejecucion de la maquina virtual de java.
	*/
	public void crearArchivosColas(){
		File temporal = new File("Temporal");
		temporal.mkdir();
		temporal.deleteOnExit();

		File queue0 = new File("Temporal\\queue0.txt");
		File queue1 = new File("Temporal\\queue1.txt");
		File queue2 = new File("Temporal\\queue2.txt");
		File queue3 = new File("Temporal\\queue3.txt");
		File queue4 = new File("Temporal\\queue4.txt");
		File queue5 = new File("Temporal\\queue5.txt");
		File queue6 = new File("Temporal\\queue6.txt");
		File queue7 = new File("Temporal\\queue7.txt");
		File queue8 = new File("Temporal\\queue8.txt");
		File queue9 = new File("Temporal\\queue9.txt");

		try{
			queue0.createNewFile();
			queue1.createNewFile();
			queue2.createNewFile();
			queue3.createNewFile();
			queue4.createNewFile();
			queue5.createNewFile();
			queue6.createNewFile();
			queue7.createNewFile();
			queue8.createNewFile();
			queue9.createNewFile();

		}catch(IOException e){
			e.printStackTrace();
		}

		queue0.deleteOnExit();
		queue1.deleteOnExit();
		queue2.deleteOnExit();
		queue3.deleteOnExit();
		queue4.deleteOnExit();
		queue5.deleteOnExit();
		queue6.deleteOnExit();
		queue7.deleteOnExit();
		queue8.deleteOnExit();
		queue9.deleteOnExit();

	}

	/**
		A partir de las 10 colas que se crearon con anterioridad en el metodo crearArchivosColas(), en este metodo se
		borra todo su contenido interno. Para volverlas a usar con el nuevo digito significativo a analizar.
	*/
	public void limpiarColas(){
		try{
			FileWriter queue0 = new FileWriter("Temporal\\queue0.txt");
			FileWriter queue1 = new FileWriter("Temporal\\queue1.txt");
			FileWriter queue2 = new FileWriter("Temporal\\queue2.txt");
			FileWriter queue3 = new FileWriter("Temporal\\queue3.txt");
			FileWriter queue4 = new FileWriter("Temporal\\queue4.txt");
			FileWriter queue5 = new FileWriter("Temporal\\queue5.txt");
			FileWriter queue6 = new FileWriter("Temporal\\queue6.txt");
			FileWriter queue7 = new FileWriter("Temporal\\queue7.txt");
			FileWriter queue8 = new FileWriter("Temporal\\queue8.txt");
			FileWriter queue9 = new FileWriter("Temporal\\queue9.txt");
			queue0.close();
			queue1.close();
			queue2.close();
			queue3.close();
			queue4.close();
			queue5.close();
			queue6.close();
			queue7.close();
			queue8.close();
			queue9.close();
		}catch(IOException e){
			e.printStackTrace();
		}
		
	}
}