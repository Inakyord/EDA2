import java.util.Scanner;
import java.io.File;

public class Principal{
	
	public static void main(String[] args){
		Scanner sc = new Scanner(System.in);

		String nombreArchivo = null;
		int opcion, opcion2, opcion3, opcion4;
		int aux=0;

		Menu menuPrincipal = new Menu(	"\n\t\tORDENAMIENTOS EXTERNOS\n"+
										"    _______________________________________________\n"+
										"    -----------------------------------------------\n"+
										"\n\n\t 1. Ingresar nombre archivo.txt\n"+
										"\n\t 2. Salir\n",2);

		Menu direccionArchivo = new Menu("\n\n\t\t     MENU ARCHIVOS\n"+
										"    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"+
										"\n\t 1. Usar archivo ingresado previamente.\n"+
										"\n\t 2. Ingresar nuevo nombre/direccion archivo.txt\n", 2);

		Menu menuOrdenamientos = new Menu(	"\n\n\t\t  MENU ORDENAMIENTOS\n"+
											"    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"+
											"\n\t 1. Polifase.\n"+
											"\n\t 2. Mezcla Equilibrada.\n"+
											"\n\t 3. Radix.\n",3);

		Menu claveOrdenamiento = new Menu(	"\n\n\t  QUE TIPO DE ORDENAMIENTO?\n"+
											"\n\t 1. Por nombre.\n"+
											"\n\t 2. Por apellidos.\n"+
											"\n\t 3. Por numero de cuenta.\n", 3);

		do{
			System.out.println("\n\n\n\n  ***************************************************");
			opcion = menuPrincipal.mostrarMenu(); //Eleccion ingresar nombre o salir.
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


			opcion3 = menuOrdenamientos.mostrarMenu(); //Eleccion poli, mezcla o radix.
			if(opcion3 == 3){ //RADIX SORT
				RadixSort ordenar = new RadixSort();
				ordenar.sort(nombreArchivo);
			}else{
				opcion4 = claveOrdenamiento.mostrarMenu(); //Elecion nombre, apellido, cuenta.
				switch(opcion3){
					case 1: //MEZCLA EQUILIBRADA
						switch(opcion4){
							case 1:
								//Mezcla equilibrada por nombre.
								break;
							case 2:
								//Mezcla equilibrada por apellidos.
								break;
							case 3:
								//Mezcla equilibrada por numero de cuenta.
								break;
						}
						break;
					case 2: //POLIFASE
						switch(opcion4){
							case 1:
								//Polifase por nombre.
								break;
							case 2:
								//Polifase por apellidos.
								break;
							case 3:
								//Polifase por numero de cuenta.
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