from PIL import Image, ImageTk
import tkinter as tk
from tkinter import ttk, filedialog
import os
import glob
import re

class ImageSplitter:
    def __init__(self, root, input_directory=None, output_directory=None):
        self.root = root
        self.root.title("Image Cropping Tool")
        
        # Setup directories
        self.setup_directories(input_directory, output_directory)
        
        # Create main frame
        self.main_frame = ttk.Frame(root)
        self.main_frame.pack(fill=tk.BOTH, expand=True)
        
        # Create left panel for controls
        self.left_panel = ttk.Frame(self.main_frame)
        self.left_panel.pack(side=tk.LEFT, fill=tk.Y, padx=10, pady=5)
        
        # Directory selection in left panel
        self.dir_frame = ttk.Frame(self.left_panel)
        self.dir_frame.pack(fill=tk.X)
        
        # Input directory selection
        ttk.Label(self.dir_frame, text="Input Dir:").pack(fill=tk.X)
        self.input_dir_var = tk.StringVar(value=self.input_dir)
        ttk.Entry(self.dir_frame, textvariable=self.input_dir_var).pack(fill=tk.X)
        ttk.Button(self.dir_frame, text="Browse", command=self.browse_input_dir).pack(fill=tk.X)
        
        # Output directory selection
        ttk.Label(self.dir_frame, text="Output Dir:").pack(fill=tk.X)
        self.output_dir_var = tk.StringVar(value=self.output_dir)
        ttk.Entry(self.dir_frame, textvariable=self.output_dir_var).pack(fill=tk.X)
        ttk.Button(self.dir_frame, text="Browse", command=self.browse_output_dir).pack(fill=tk.X)
        
        # File name entry
        ttk.Label(self.left_panel, text="Output Name:").pack(fill=tk.X, pady=(10,2))
        self.output_name_var = tk.StringVar()
        ttk.Entry(self.left_panel, textvariable=self.output_name_var).pack(fill=tk.X)
        
        # Control buttons
        ttk.Button(self.left_panel, text="Add Horizontal Line",
                  command=self.add_horizontal_line).pack(fill=tk.X, pady=(10,2))
        ttk.Button(self.left_panel, text="Add Vertical Line",
                  command=self.add_vertical_line).pack(fill=tk.X, pady=2)
        ttk.Button(self.left_panel, text="Next Image",
                  command=self.next_image).pack(fill=tk.X, pady=2)
        ttk.Button(self.left_panel, text="Previous Image",
                  command=self.prev_image).pack(fill=tk.X, pady=2)
        ttk.Button(self.left_panel, text="Crop Image",
                  command=self.split_image).pack(fill=tk.X, pady=2)
        
        # Create canvas frame
        self.canvas_frame = ttk.Frame(self.main_frame)
        self.canvas_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        
        # Create canvas
        self.canvas = tk.Canvas(self.canvas_frame)
        self.canvas.pack(fill=tk.BOTH, expand=True)
        
        # Initialize cutting lines
        self.h_lines = []  # horizontal lines
        self.v_lines = []  # vertical lines
        self.active_line = None
        
        # Bind events
        self.canvas.bind('<Button-1>', self.on_click)
        self.canvas.bind('<B1-Motion>', self.on_drag)
        self.canvas.bind('<ButtonRelease-1>', self.on_release)
        
        # Load images if directory exists
        self.load_image_list()
        
    def setup_directories(self, input_directory, output_directory):
        """Setup input and output directories"""
        self.input_dir = input_directory or ""
        self.output_dir = output_directory or ""
        
    def browse_input_dir(self):
        """Browse for input directory"""
        directory = filedialog.askdirectory()
        if directory:
            self.input_dir = directory
            self.input_dir_var.set(directory)
            self.load_image_list()
            
    def browse_output_dir(self):
        """Browse for output directory"""
        directory = filedialog.askdirectory()
        if directory:
            self.output_dir = directory
            self.output_dir_var.set(directory)
            
    def load_image_list(self):
        """Load and sort image list"""
        if not os.path.exists(self.input_dir):
            return
            
        # Get all PNG files
        self.image_files = glob.glob(os.path.join(self.input_dir, "*.png"))
        
        # Sort files naturally (e.g., 1, 2, 10 instead of 1, 10, 2)
        self.image_files.sort(key=lambda x: [int(c) if c.isdigit() else c.lower() for c in re.split('([0-9]+)', x)])
        
        self.current_image_index = 0
        
        # Load first image if available
        if self.image_files:
            self.load_image(self.image_files[0])
            # Set default output name
            base_name = os.path.splitext(os.path.basename(self.image_files[0]))[0]
            self.output_name_var.set(f"{base_name}")
    
    def load_image(self, image_path):
        """Load a new image"""
        if not os.path.exists(image_path):
            return
            
        self.original_image = Image.open(image_path)
        self.current_image_path = image_path
        
        # Update output name
        base_name = os.path.splitext(os.path.basename(image_path))[0]
        self.output_name_var.set(f"{base_name}")
        
        self.update_display()
    
    def update_display(self):
        """Update the display"""
        self.display_width = self.original_image.width
        self.display_height = self.original_image.height
        self.display_image = self.original_image
        
        self.photo = ImageTk.PhotoImage(self.display_image)
        self.canvas.delete("all")
        self.canvas.config(width=self.display_width, height=self.display_height)
        self.canvas.create_image(0, 0, image=self.photo, anchor=tk.NW)
        
        # Add default lines after loading image
        self.add_default_lines()
    
    def redraw_lines(self):
        """Redraw all lines"""
        old_h_lines = self.h_lines.copy()
        old_v_lines = self.v_lines.copy()
        
        self.h_lines = []
        self.v_lines = []
        
        for line in old_h_lines:
            coords = self.canvas.coords(line)
            if coords:
                y = coords[1]
                self.add_horizontal_line(y)
        
        for line in old_v_lines:
            coords = self.canvas.coords(line)
            if coords:
                x = coords[0]
                self.add_vertical_line(x)
    
    def next_image(self):
        """Load next image in the directory"""
        if not self.image_files:
            return
        self.current_image_index = (self.current_image_index + 1) % len(self.image_files)
        self.load_image(self.image_files[self.current_image_index])
    
    def prev_image(self):
        """Load previous image in the directory"""
        if not self.image_files:
            return
        self.current_image_index = (self.current_image_index - 1) % len(self.image_files)
        self.load_image(self.image_files[self.current_image_index])
    
    def add_horizontal_line(self, y=None):
        """Add horizontal cutting line"""
        if y is None:
            y = self.display_height // 2
        line = self.canvas.create_line(0, y, self.display_width, y,
                                     fill='red', width=2)
        self.h_lines.append(line)
    
    def add_vertical_line(self, x=None):
        """Add vertical cutting line"""
        if x is None:
            x = self.display_width // 2
        line = self.canvas.create_line(x, 0, x, self.display_height,
                                     fill='blue', width=2)
        self.v_lines.append(line)
    
    def on_click(self, event):
        """Handle click events"""
        x, y = event.x, event.y
        for line in self.h_lines + self.v_lines:
            coords = self.canvas.coords(line)
            if self.is_near_line(x, y, coords):
                self.active_line = line
                break
    
    def is_near_line(self, x, y, coords):
        """Check if click position is near a line"""
        if len(coords) == 4:  # line coordinates: x1,y1,x2,y2
            x1, y1, x2, y2 = coords
            if abs(x1 - x2) < 2:  # vertical line
                return abs(x - x1) < 5
            if abs(y1 - y2) < 2:  # horizontal line
                return abs(y - y1) < 5
        return False
    
    def on_drag(self, event):
        """Handle drag events"""
        if self.active_line:
            if self.active_line in self.h_lines:
                self.canvas.coords(self.active_line,
                                 0, event.y,
                                 self.display_width, event.y)
            else:
                self.canvas.coords(self.active_line,
                                 event.x, 0,
                                 event.x, self.display_height)
    
    def on_release(self, event):
        """Handle release events"""
        self.active_line = None
    
    def split_image(self):
        """Crop the image to keep only the center part"""
        if len(self.h_lines) < 2 or len(self.v_lines) < 2:
            print("請至少添加2條水平線和2條垂直線")
            return
            
        if not self.output_dir:
            print("請選擇輸出目錄")
            return
            
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        
        # Get line positions
        h_positions = sorted([self.canvas.coords(line)[1] for line in self.h_lines])
        v_positions = sorted([self.canvas.coords(line)[0] for line in self.v_lines])
        
        # Get middle section coordinates
        middle_h = len(h_positions) // 2
        middle_v = len(v_positions) // 2
        
        # Calculate crop coordinates for middle section
        left = int(v_positions[middle_v - 1])
        upper = int(h_positions[middle_h - 1])
        right = int(v_positions[middle_v])
        lower = int(h_positions[middle_h])
        
        # Crop image
        section = self.original_image.crop((left, upper, right, lower))
        
        # Get output filename
        output_name = self.output_name_var.get() or "cropped_image"
        if not output_name.endswith('.png'):
            output_name += '.png'
            
        # Save cropped image
        output_path = os.path.join(self.output_dir, output_name)
        section.save(output_path)
        
        # Save cutting line parameters
        txt_output_name = os.path.splitext(output_name)[0] + '.txt'
        txt_output_path = os.path.join(self.output_dir, txt_output_name)
        
        with open(txt_output_path, 'w') as f:
            f.write("Horizontal Lines (y-coordinates):\n")
            for pos in h_positions:
                f.write(f"{int(pos)}\n")
            f.write("\nVertical Lines (x-coordinates):\n")
            for pos in v_positions:
                f.write(f"{int(pos)}\n")
            f.write("\nCrop Coordinates:\n")
            f.write(f"Left: {left}\n")
            f.write(f"Upper: {upper}\n")
            f.write(f"Right: {right}\n")
            f.write(f"Lower: {lower}\n")
        
        print(f"已保存裁剪後的圖片到: {output_path}")
        print(f"已保存切割線參數到: {txt_output_path}")

    def add_default_lines(self):
        """Add default cutting lines"""
        if not hasattr(self, 'display_width') or not hasattr(self, 'display_height'):
            return
            
        # Clear existing lines
        self.h_lines = []
        self.v_lines = []
        
        # Add two horizontal lines at 1/3 and 2/3 height
        h1 = self.display_height // 3
        h2 = (self.display_height * 2) // 3
        self.add_horizontal_line(h1)
        self.add_horizontal_line(h2)
        
        # Add two vertical lines at 1/3 and 2/3 width
        v1 = self.display_width // 3
        v2 = (self.display_width * 2) // 3
        self.add_vertical_line(v1)
        self.add_vertical_line(v2)

def run_image_splitter(input_dir=None, output_dir=None):
    """Run the image splitter application"""
    root = tk.Tk()
    app = ImageSplitter(root, input_dir, output_dir)
    root.mainloop()

if __name__ == "__main__":
    run_image_splitter() 