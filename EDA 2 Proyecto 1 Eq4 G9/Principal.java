import java.util.Scanner;
import java.io.File;
import mezclaequilibrada.codigos.*;
import polifase.*;

/**
 * Clase principal del programa para el proyecto 1 de EDA 2.
 * @author David Calderon Jimenez
 * @author Humberto Ignacio Hernandez Olvera
 * @author I\u00f1aky Ordiales Caballero
 * @version 1.7
 * @since 22/Octubre/2020
 */
public class Principal{

	public static void main(String[] args){
		Scanner sc = new Scanner(System.in);

		String nombreArchivo = null;
		int opcion, opcion2, opcion3, opcion4;
		int aux=0;


		/* Esta instancia de la clase Menu es el menu inicial donde se elige si ingresar archivo para ordenar
		   o salir del programa.*/
		Menu menuPrincipal = new Menu(	"\n\t\tORDENAMIENTOS EXTERNOS\n"+
										"    _______________________________________________\n"+
										"    -----------------------------------------------\n"+
										"\n\n\t 1. Ingresar nombre archivo.txt\n"+
										"\n\t 2. Salir\n",2);

		/* Esta instancia de la clase Menu es el menu donde en caso de ya haber ingresado un archivo previamente,
		   te da la opcion de volver a usar ese archivo o de ingresar uno nuevo.*/
		Menu direccionArchivo = new Menu("\n\n\t\t     MENU ARCHIVOS\n"+
										"    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"+
										"\n\t 1. Usar archivo ingresado previamente.\n"+
										"\n\t 2. Ingresar nuevo nombre/direccion archivo.txt\n", 2);

		/* Esta instancia de la clase Menu te da las opciones de los diferentes algoritmos de ordenamiento
		   externo: Polifase, Mezcla o Radix.*/
		Menu menuOrdenamientos = new Menu(	"\n\n\t\t  MENU ORDENAMIENTOS\n"+
											"    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"+
											"\n\t 1. Polifase.\n"+
											"\n\t 2. Mezcla Equilibrada.\n"+
											"\n\t 3. Radix.\n",3);

		/* Esta instancia de la clase Menu te da la opcion de ordenar por nombre, apellidos o numero de cuenta*/
		Menu claveOrdenamiento = new Menu(	"\n\n\t  QUE TIPO DE ORDENAMIENTO?\n"+
											"\n\t 1. Por nombre.\n"+
											"\n\t 2. Por apellidos.\n", 2);

		// A continuacion se empiezan a mostrar los menus para el programa.
		do{
			System.out.println("\n\n\n\n  ***************************************************");
			opcion = menuPrincipal.mostrarMenu(); //Eleccion ingresar archivo o salir.
			switch(opcion){
				case 1:
					if(nombreArchivo!=null){
						// POSIBILIDAD DE USAR ANTIGUO ARCHIVO.
						opcion2 = direccionArchivo.mostrarMenu(); //Eleccion usar archivo ya ingresado, ingresar nuevo.
					}else
						opcion2 = 2;

					if(opcion2 == 2){
						// INGRESAR NUEVO NOMBRE ARCHIVO
						System.out.print("\n\t Ingrese nombre de Archivo con su extension (.txt al final): ");
						nombreArchivo = sc.nextLine();
						while (verificar(nombreArchivo)!=true){
							aux = 0;
							do{
								System.out.print("\n\t Desea volver a ingresar nombre de archivo?:\n\t  1. Si quiero.\n\t  2. No, Salir.\n\t Opcion: ");
								opcion2 = sc.nextInt();
								switch (opcion2){
									case 1:
										System.out.print("\n\t Ingrese nombre archivo (.txt al final): ");
										sc.nextLine();
										nombreArchivo = sc.nextLine();
										break;
									case 2:
										sc.nextLine();
										aux = 2;
										break;
									default:
										System.out.println("\n\t Opcion no valida, vuelvalo a intentar.");
										break;
								}
							}while(opcion2 != 2);
							if (aux == 2){
								opcion = 2;
								break;
							}
						}//Cierra while de verificacion
					}
					break; // Acaba caso 1.
				case 2:
					break;
			} // Cierra switch sobre el nombre del archivo.
			if (opcion==2){
				//Se sale del programa si en algun momento se escogió la opcion de salir.
				break;
			}


			//A partir de aqui empiezan los ordenamientos. |||||


			opcion3 = menuOrdenamientos.mostrarMenu(); //Eleccion polifase, mezcla o radix.
			if(opcion3 == 3){ //RADIX SORT
				RadixSort ordenar = new RadixSort();
				ordenar.sort(nombreArchivo);
			}else{
				opcion4 = claveOrdenamiento.mostrarMenu(); //Elecion nombre, apellido, cuenta.
				switch(opcion3){
					case 2: //MEZCLA EQUILIBRADA solo por nombre o apellido
						switch(opcion4){
							case 1:	// Mezcla equilibrada por nombre.
								MezclaEquilibrada meq = new MezclaEquilibrada();
								meq.ordenar(nombreArchivo, false);
								break;
							case 2:	//Mezcla equilibrada por apellidos.
								meq = new MezclaEquilibrada();
								meq.ordenar(nombreArchivo, true);
								break;
						}
						break;
					case 1: //POLIFASE
						switch(opcion4){
							case 1://Polifase por nombres
								Polifase poli = new Polifase();
								poli.ordenar(nombreArchivo, false);
								break;
							case 2://Polifase por apellidos
								poli = new Polifase();
								poli.ordenar(nombreArchivo, true);
								break;
						}
						break;
				}
			}
			System.out.println("\n\n\t EL ARCHIVO \""+nombreArchivo+"\" SE ORDENO CON EXITO\n\n\t Presione Enter para volver al menu principal.");
			sc.nextLine();
		}while (true);

	System.out.println("\n\n\t *** *** ADIOS ...\n");

	} //Cierra main





	/**
		El metodo verificar recibe el nombre/direccion de un archivo. Dentro del metodo se verifica si se encuentra con
		exito el archivo o no. En caso de encontrarlo regresa true, si no lo encuentra regresa false.
		@param nombreArchivo El nombre del archivo que se ordenara.
		@return true si el archivo existe, false si el archivo no existe.
	*/
	public static boolean verificar(String nombreArchivo){
		File archivo = new File(nombreArchivo);
		if (archivo.exists()){
			System.out.println("\n\t El archivo se encontro con exito.");
			return true;
		}else{
			System.out.println("\n\t No se encuentra el archivo !!!");
			return false;
		}

	}
}
