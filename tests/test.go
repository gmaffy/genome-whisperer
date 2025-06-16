package main

import (
	"fmt"
	"os"
	"os/exec"
)

func main() {
	fmt.Println("Hello World")
	err := os.MkdirAll("masurca", 0755)
	if err != nil {
		return
	}
	os.Chdir("masurca")
	cmd2 := exec.Command("touch", "test.txt")
	err = cmd2.Run()
	if err != nil {
		fmt.Println("CMD error:", err)
		return
	}
	cmd3 := exec.Command("ls")
	err = cmd3.Run()
	if err != nil {
		fmt.Println("CMD error:", err)
	}
	fmt.Println("Done")

}
